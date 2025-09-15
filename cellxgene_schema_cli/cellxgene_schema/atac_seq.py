import contextlib
import gzip
import logging
import os
import shutil
import subprocess as sp
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import anndata as ad
import h5py
import ibis
import pyarrow as pa
import pyarrow.csv
import pyarrow.dataset
import pysam

from .ontology_parser import ONTOLOGY_PARSER
from .utils import GB, is_ontological_descendant_of

logger = logging.getLogger(__name__)

# TODO: these chromosome tables should be calculated from the fasta file?
# location of fasta https://www.gencodegenes.org/human/release_44.html and file name GRCh38.primary_assembly.genome.fa
human_chromosome_by_length = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569,
    "GL000009.2": 201709,
    "GL000194.1": 191469,
    "GL000195.1": 182896,
    "GL000205.2": 185591,
    "GL000213.1": 164239,
    "GL000216.2": 176608,
    "GL000218.1": 161147,
    "GL000219.1": 179198,
    "GL000220.1": 161802,
    "GL000225.1": 211173,
    "KI270442.1": 392061,
    "KI270711.1": 42210,
    "KI270713.1": 40745,
    "KI270721.1": 100316,
    "KI270726.1": 43739,
    "KI270727.1": 448248,
    "KI270728.1": 1872759,
    "KI270731.1": 150754,
    "KI270733.1": 179772,
    "KI270734.1": 165050,
    "KI270744.1": 168472,
    "KI270750.1": 148850,
}
mouse_chromosome_by_length = {
    "chr1": 195154279,
    "chr2": 181755017,
    "chr3": 159745316,
    "chr4": 156860686,
    "chr5": 151758149,
    "chr6": 149588044,
    "chr7": 144995196,
    "chr8": 130127694,
    "chr9": 124359700,
    "chr10": 130530862,
    "chr11": 121973369,
    "chr12": 120092757,
    "chr13": 120883175,
    "chr14": 125139656,
    "chr15": 104073951,
    "chr16": 98008968,
    "chr17": 95294699,
    "chr18": 90720763,
    "chr19": 61420004,
    "chrX": 169476592,
    "chrY": 91455967,
    "chrM": 16299,
    "GL456210.1": 169725,
    "GL456211.1": 241735,
    "GL456212.1": 153618,
    "GL456219.1": 175968,
    "GL456221.1": 206961,
    "GL456239.1": 40056,
    "GL456354.1": 195993,
    "GL456372.1": 28664,
    "GL456381.1": 25871,
    "GL456385.1": 35240,
    "JH584295.1": 1976,
    "JH584296.1": 199368,
    "JH584297.1": 205776,
    "JH584298.1": 184189,
    "JH584299.1": 953012,
    "JH584303.1": 158099,
    "JH584304.1": 114452,
}
organism_ontology_term_id_by_chromosome_length_table = {
    "NCBITaxon:9606": human_chromosome_by_length,
    "NCBITaxon:10090": mouse_chromosome_by_length,
}
column_ordering = ["chromosome", "start coordinate", "stop coordinate", "barcode", "read support"]
schema = pa.schema(
    [
        pa.field("chromosome", pa.string()),
        pa.field("start coordinate", pa.int64()),
        pa.field("stop coordinate", pa.int64()),
        pa.field("barcode", pa.string()),
        pa.field("read support", pa.int64()),
    ]
)


def log_calls(func):
    def wrapper(*args, **kwargs):
        logging.info("starting: %s", func.__name__)
        return func(*args, **kwargs)

    return wrapper


def count_lines_in_compressed_file(file_path: str) -> int:
    """
    Count lines in a compressed file (gzip or bgzip).

    :param file_path: Path to the compressed file
    :return: Number of lines in the file
    """
    line_count = 0
    with gzip.open(file_path, "rt") as f:
        for _ in f:
            line_count += 1
    return line_count


def is_atac(x: str) -> str:
    if is_ontological_descendant_of(ONTOLOGY_PARSER, x, "EFO:0010891"):
        if is_ontological_descendant_of(ONTOLOGY_PARSER, x, "EFO:0008913"):
            return "p"  # paired
        else:
            return "u"  # unpaired
    else:
        return "n"  # not atac seq


def check_anndata_requires_fragment(anndata_file: str) -> bool:
    """
    Check if an anndata file requires a fragment file to be valid. The anndata file requires a fragment file if the
    assay_ontology_term_id are all descendants of EFO:0010891 and EFO:0008913. If the assay_ontology_term_id are all
    descendants of EFO:0010891 and not EFO:0008913, the anndata file does not require a fragment file. If the
    assay_ontology_term_id are not all descendants of EFO:0010891, an error is raised and a fragment file is not allowed
    because the anndata is not a valid ATAC-Seq file.

    :param anndata_file: The anndata file to validate.
    :return:
    """
    with h5py.File(anndata_file) as f:
        assay_ontology_term_ids = ad.io.read_elem(f["obs"])["assay_ontology_term_id"]
    df = assay_ontology_term_ids.map(is_atac)

    if (df == "p").all():
        return False
    elif (df == "u").all():
        return True
    if (df == "n").any():
        raise ValueError("Anndata.obs.assay_ontology_term_id are not all descendants of EFO:0010891.")
    else:
        raise ValueError("Anndata.obs.assay_ontology_term_id has mixed paired and unpaired assay terms.")


def process_fragment(
    fragment_file: str,
    anndata_file: str,
    generate_index: bool = False,
    output_file: Optional[str] = None,
) -> list[str]:
    """
    Validate the fragment against the anndata file and generate the index if the fragment is valid.

    :param str fragment_file: The fragment file to process
    :param str anndata_file: The anndata file to validate against
    :param bool generate_index: Whether to generate the index for the fragment
    :param override_write_algorithm: Override the write algorithm used to write the bgzip file. Options are "pysam"
    and "cli"
    :param output_file: The output file to write the bgzip file to. If not provided, the output file will be the same

    """

    with tempfile.TemporaryDirectory() as tempdir:
        # quick checks
        errors = validate_anndata(anndata_file)
        if errors:
            return errors

        # convert the fragment to a parquet file for faster processing
        try:
            parquet_file = convert_to_parquet(fragment_file, tempdir)
        except Exception as e:
            msg = "Error Parsing the fragment file. Check that columns match schema definition. Error: " + str(e)
            logger.exception(msg)
            return [msg]

        organism_ontology_term_id = None
        with h5py.File(anndata_file) as f:
            organism_ontology_term_id = ad.io.read_elem(f["uns"]).get("organism_ontology_term_id")

        # slow checks
        errors = validate_anndata_with_fragment(parquet_file, anndata_file, organism_ontology_term_id)
        if errors:
            return errors
        else:
            logger.info("Fragment and Anndata file are valid")

    # generate the index
    if generate_index:
        logger.info(f"Sorting fragment and generating index for {fragment_file}")
        index_fragment(fragment_file, output_file)

    logger.debug("cleaning up")
    return []


def convert_to_parquet(fragment_file: str, tempdir: str) -> str:
    """
    Convert the fragment file to a parquet dataset for faster processing.

    :param fragment_file: A gzipped compressed fragment file
    :param tempdir: The temporary directory to write the parquet file to. Name of the written file is derived from
    the input.
    """
    logger.info(f"Converting {fragment_file} to parquet")
    parquet_file_path = Path(tempdir) / Path(fragment_file).name.replace(".gz", ".parquet").replace(".bgz", ".parquet")
    pa.dataset.write_dataset(
        data=pa.csv.open_csv(
            pa.input_stream(fragment_file, compression="gzip", buffer_size=GB),
            read_options=pa.csv.ReadOptions(column_names=schema.names),
            parse_options=pa.csv.ParseOptions(delimiter="\t"),
            convert_options=pa.csv.ConvertOptions(column_types=schema),
        ),
        base_dir=parquet_file_path,
        format="parquet",
        # Using hive partitioning for best compatibility with dask to_parquet and read_parquet functions
        partitioning=pa.dataset.partitioning(pa.schema([pa.field("chromosome", pa.string())]), flavor="hive"),
    )

    return str(parquet_file_path)


def report_errors(header: str, errors: list[str]) -> list[str]:
    if any(errors):
        errors = [e for e in errors if e is not None]
        errors = [header] + errors
        logger.error("\n\t".join(errors))
        return errors
    else:
        return []


def validate_anndata(anndata_file: str) -> list[str]:
    errors = [validate_anndata_organism_ontology_term_id(anndata_file), validate_anndata_is_primary_data(anndata_file)]
    return report_errors("Errors found in Anndata file. Skipping fragment validation.", errors)


def validate_anndata_with_fragment(parquet_file: str, anndata_file: str, organism_ontology_term_id: str) -> list[str]:
    errors = [
        validate_fragment_start_coordinate_greater_than_0(parquet_file),
        validate_fragment_barcode_in_adata_index(parquet_file, anndata_file),
        validate_fragment_stop_coordinate_within_chromosome(parquet_file, organism_ontology_term_id),
        validate_fragment_stop_greater_than_start_coordinate(parquet_file),
        validate_fragment_read_support(parquet_file),
        validate_fragment_no_duplicate_rows(parquet_file),
    ]
    return report_errors("Errors found in Fragment and/or Anndata file", errors)


@log_calls
def validate_fragment_no_duplicate_rows(parquet_file: str) -> Optional[str]:
    t = ibis.read_parquet(f"{parquet_file}/**", hive_partitioning=True)
    rows_per_chromosome = t["chromosome"].value_counts().execute()
    msg = ""
    # Checking number of unique rows per chromosome is more memory efficient than checking all rows at once
    for chromosome, count in rows_per_chromosome.itertuples(index=False):
        n_unique = t.filter(t["chromosome"] == chromosome).distinct().count().execute()
        if n_unique != count:
            if not msg:
                msg = "Fragment file has duplicate rows.\n"
            msg += f"Chromosome {chromosome} has {count} rows but only {n_unique} are unique\n"
    if msg:
        return msg.strip()  # remove trailing newline


@log_calls
def validate_fragment_start_coordinate_greater_than_0(parquet_file: str) -> Optional[str]:
    df = ibis.read_parquet(f"{parquet_file}/**", hive_partitioning=True)
    if not (df["start coordinate"] > 0).all().execute():
        return "Start coordinate must be greater than 0."


@log_calls
def validate_fragment_barcode_in_adata_index(parquet_file: str, anndata_file: str) -> Optional[str]:
    df = ibis.read_parquet(f"{parquet_file}/**", hive_partitioning=True)
    barcode = set(df.select("barcode").distinct().execute()["barcode"])
    with h5py.File(anndata_file) as f:
        obs = ad.io.read_elem(f["obs"])
    if set(obs.index) != barcode:
        return "Barcodes don't match anndata.obs.index"


@log_calls
def validate_fragment_stop_greater_than_start_coordinate(parquet_file: str) -> Optional[str]:
    df = ibis.read_parquet(f"{parquet_file}/**", hive_partitioning=True)
    if not (df["stop coordinate"] > df["start coordinate"]).all().execute():
        return "Stop coordinate must be greater than start coordinate."


@log_calls
def validate_fragment_stop_coordinate_within_chromosome(
    parquet_file: str, organism_ontology_term_id: str
) -> Optional[str]:
    chromosome_length_table = organism_ontology_term_id_by_chromosome_length_table.get(organism_ontology_term_id)
    t = ibis.read_parquet(f"{parquet_file}/**", hive_partitioning=True)
    df = t.group_by("chromosome").aggregate(max_stop_coordinate=t["stop coordinate"].max()).execute()

    mismatched_chromosomes = set(df["chromosome"].unique()) - chromosome_length_table.keys()
    if mismatched_chromosomes:
        return f"Chromosomes in the fragment do not match the organism({organism_ontology_term_id}).\n" + "\t\n".join(
            mismatched_chromosomes
        )

    df["chromosome_length"] = df["chromosome"].map(chromosome_length_table)
    if not (df["max_stop_coordinate"] <= df["chromosome_length"]).all():
        return "Stop coordinate must be less than the chromosome length."


@log_calls
def validate_fragment_read_support(parquet_file: str) -> Optional[str]:
    t = ibis.read_parquet(f"{parquet_file}/**", hive_partitioning=True)
    if (t["read support"] <= 0).any().execute():
        return "Read support must be greater than 0."


@log_calls
def validate_anndata_is_primary_data(anndata_file: str) -> Optional[str]:
    with h5py.File(anndata_file) as f:
        is_primary_data = ad.io.read_elem(f["obs"])["is_primary_data"]
    if not is_primary_data.all():
        return "Anndata.obs.is_primary_data must all be True."


@log_calls
def validate_anndata_organism_ontology_term_id(anndata_file: str) -> Optional[str]:
    organism_ontology_term_id = None
    with h5py.File(anndata_file) as f:
        organism_ontology_term_id = ad.io.read_elem(f["uns"]).get("organism_ontology_term_id")
    allowed_terms = [*organism_ontology_term_id_by_chromosome_length_table.keys()]
    if organism_ontology_term_id not in allowed_terms:
        return f"Anndata.obs.organism_ontology_term_id must be one of {allowed_terms}. Got {organism_ontology_term_id}."


@log_calls
def detect_chromosomes(parquet_file: str) -> list[str]:
    t = ibis.read_parquet(f"{parquet_file}/**", hive_partitioning=True)
    chromosomes = list(t.select(["chromosome"]).distinct().execute()["chromosome"])
    chromosomes.sort()  # sort chromosomes to ensure consistent order
    return chromosomes


def get_output_file(fragment_file: str, output_file: Optional[str]) -> str:
    if not output_file:
        bgzip_output_file = fragment_file.replace(".gz", ".bgz")
    elif not output_file.endswith(".bgz"):
        bgzip_output_file = output_file + ".bgz"
    else:
        bgzip_output_file = output_file
    return bgzip_output_file


def index_fragment(
    fragment_file: str,
    output_file: Optional[str] = None,
):
    # sort the fragment by chromosome, start coordinate, and stop coordinate, then compress it with bgzip
    bgzip_output_file = get_output_file(fragment_file, output_file)
    bgzip_output_path = Path(bgzip_output_file)
    bgzip_output_path.unlink(missing_ok=True)

    prepare_fragment(fragment_file, bgzip_output_file)
    logger.info(f"Fragment sorted and compressed: {bgzip_output_file}")
    pysam.tabix_index(bgzip_output_file, preset="bed", force=True)
    tabix_output_file = bgzip_output_file + ".tbi"
    logger.info(f"Index file generated: {tabix_output_file}")

    # Validate line count between original and output files
    logger.info("Validating line count between original and output files")

    # Count lines in parallel to reduce processing time for large files
    with ThreadPoolExecutor(max_workers=2) as executor:
        original_future = executor.submit(count_lines_in_compressed_file, fragment_file)
        output_future = executor.submit(count_lines_in_compressed_file, bgzip_output_file)

        original_line_count = original_future.result()
        output_line_count = output_future.result()

    if original_line_count != output_line_count:
        error_msg = f"Line count validation failed: original file has {original_line_count} lines, output file has {output_line_count} lines"
        logger.error(error_msg)
        raise ValueError(error_msg)

    logger.info(f"Line count validation passed: {original_line_count} lines in both original and output files")


SORT_MEMORY_PERCENTAGE = 50


def _calculate_sort_memory(num_cores: int, sort_memory_percent: int) -> int:
    """
    Calculate the memory percentage to allocate per sort thread.

    :param num_cores: Number of CPU cores available.
    :param sort_memory_percent: Total percentage of system memory to allocate for sorting.
    :return: Percentage of memory to allocate per sort thread.
    """
    effective_memory_pct = min(sort_memory_percent, 50) if num_cores > 1 else sort_memory_percent
    return max(effective_memory_pct // num_cores, 1)  # ensure at least 1% memory per core


def _ensure_tools(tools: Iterable[Tuple[str, str]]) -> None:
    """Ensure required CLI tools are present in PATH."""
    missing = [name for name, _ in tools if shutil.which(name) is None]
    if missing:
        msgs = dict(tools)
        raise RuntimeError(" | ".join(f"Missing '{name}' in PATH; {msgs[name]}." for name in missing))


def _default_cores() -> int:
    """Get CPU count, respecting container limits."""
    # Try to read cgroup CPU limit (Docker/Kubernetes) - most reliable method
    try:
        # cgroup v1
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as f:
            quota = int(f.read().strip())
        with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as f:
            period = int(f.read().strip())
        if quota > 0 and period > 0:
            return max(1, quota // period)
    except (FileNotFoundError, ValueError, PermissionError):
        pass

    try:
        # cgroup v2
        with open("/sys/fs/cgroup/cpu.max") as f:
            content = f.read().strip()
        if content != "max":
            parts = content.split()
            if len(parts) == 2:
                quota, period = int(parts[0]), int(parts[1])
                return max(1, quota // period)
    except (FileNotFoundError, ValueError, PermissionError):
        pass

    # Fallback to system CPU count
    return max(1, (os.cpu_count() or 1))


def _deterministic_env() -> dict:
    # GNU coreutils behaviors (sort/uniq) are locale-sensitive; pin for determinism.
    return {**os.environ, "LC_ALL": "C"}


def _sort_command(num_cores: int, sort_mem_pct: int) -> List[str]:
    """Tab-delimited sort by (chrom, start, end, strand), numeric where needed, parallel, compressed temps via pigz."""
    return [
        "sort",
        "--parallel",
        str(num_cores),
        "-t",
        "\t",
        "-k1,1",
        "-k2,2n",
        "-k3,3n",
        "-k4,4",
        "-S",
        f"{sort_mem_pct}%",
        "--compress-program",
        "pigz",
    ]


def _bgzip_command(num_cores: int) -> List[str]:
    # Prefer long form; widely supported on recent htslib/bgzip.
    return ["bgzip", f"--threads={num_cores}", "-c"]


def _pigz_decompress_command(path: Path) -> List[str]:
    return ["pigz", "-dc", str(path)]


def _pipeline_run(
    stages: Sequence[Tuple[str, List[str]]],
    output_file: Path,
    *,
    env: Optional[dict] = None,
) -> None:
    """
    Run a left-to-right streaming pipeline and write the last stage's stdout to output_file.
    Each stage is (name, argv). Captures stderr per stage for error reporting.
    """
    if not stages:
        raise ValueError("pipeline must have at least one stage")

    procs: List[sp.Popen] = []

    try:
        # open sink
        with open(output_file, "wb") as out_f:
            # launch first
            name0, cmd0 = stages[0]
            logger.info(f"Starting pipeline stage: {name0}")
            p0 = sp.Popen(cmd0, stdout=sp.PIPE, stderr=sp.PIPE, env=env)
            procs.append(p0)

            # chain middle (if any)
            prev_stdout = p0.stdout
            for stage_name, stage_cmd in stages[1:-1]:
                logger.info(f"Starting pipeline stage: {stage_name}")
                p = sp.Popen(stage_cmd, stdin=prev_stdout, stdout=sp.PIPE, stderr=sp.PIPE, env=env)
                procs.append(p)
                if prev_stdout is not None:
                    prev_stdout.close()  # allow upstream SIGPIPE if we fail downstream
                prev_stdout = p.stdout

            # last -> file
            last_name, last_cmd = stages[-1]
            logger.info(f"Starting pipeline stage: {last_name}")
            plast = sp.Popen(last_cmd, stdin=prev_stdout, stdout=out_f, stderr=sp.PIPE, env=env)
            procs.append(plast)
            if prev_stdout is not None:
                prev_stdout.close()

            plast.wait()

            # collect return codes (downstream-first)
            returncodes = []
            for proc in reversed(procs):
                rc = proc.returncode if proc is plast else proc.wait()
                returncodes.append(rc)

            # error aggregation
            if any(rc != 0 for rc in returncodes):
                errs = []

                for (stage_name, _), proc in zip(stages, procs):
                    if proc.returncode != 0:
                        try:
                            err = proc.stderr.read().decode(errors="replace").strip() if proc.stderr else ""
                        except Exception:
                            err = ""

                        errs.append(f"{stage_name} failed (rc={proc.returncode}): {err}")

                # Provide AWS Batch specific guidance
                error_msg = " | ".join(errs)

                raise RuntimeError(error_msg)

    finally:
        # harden against leaked children in edge cases
        for p in procs:
            if p.poll() is None:
                with contextlib.suppress(Exception):
                    p.terminate()


def prepare_fragment(fragment_file: str, bgzip_output_file: str) -> None:
    """
    Decompress -> sort -> bgzip.
    Sorts by (chrom, start, end, strand). Output is bgzipped; chromosomes may appear in non-input order,
    but are internally sorted within chromosome.
    """
    _ensure_tools(
        [
            ("sort", "required to sort the fragment file"),
            ("bgzip", "required to bgzip-compress the result"),
            ("pigz", "required for parallel gzip (de)compression"),
        ]
    )

    in_path = Path(fragment_file)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path}")

    out_path = Path(bgzip_output_file)
    ncores = _default_cores()
    sort_mem_pct = _calculate_sort_memory(ncores, SORT_MEMORY_PERCENTAGE)
    env = _deterministic_env()

    stages = [
        ("pigz -dc", _pigz_decompress_command(in_path)),
        ("sort", _sort_command(ncores, sort_mem_pct)),
        ("bgzip", _bgzip_command(ncores)),
    ]

    _pipeline_run(stages, out_path, env=env)
    logger.info(f"bgzip compression completed successfully for {out_path}")


def deduplicate_fragment_rows(
    fragment_file_name: str,
    output_file_name: Optional[str] = None,
    sort_memory_percent: int = SORT_MEMORY_PERCENTAGE,
) -> str:
    """
    Decompress -> sort -> uniq -> bgzip.
    Removes duplicate lines after sorting by (chrom, start, end, strand).

    :param fragment_file_name: Path to the input gzipped fragment file.
    :param output_file_name: Path for the output deduplicated bgzipped file. If None, a default name is generated.
    :param sort_memory_percent: Percentage of system memory to allocate per sort thread.
    :return: Path to the deduplicated bgzipped output file.
    :raises RuntimeError: If the operation fails due to resource constraints.
    """
    _ensure_tools(
        [
            ("sort", "required to sort the fragment file"),
            ("bgzip", "required to bgzip-compress the result"),
            ("pigz", "required for parallel gzip (de)compression"),
            ("uniq", "required to remove duplicate rows"),
        ]
    )

    in_path = Path(fragment_file_name)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path}")

    if output_file_name is None:
        # <input minus .tsv[.gz]>_dedup.tsv.bgz
        stem_once = in_path.stem  # drops final suffix (e.g., ".gz")
        stem_twice = Path(stem_once).stem  # drops ".tsv" if present
        out_path = in_path.with_name(f"{stem_twice}_dedup.tsv.bgz")
    else:
        out_path = Path(output_file_name)

    ncores = _default_cores()
    sort_mem_pct = _calculate_sort_memory(ncores, sort_memory_percent)  # Output is smaller

    stages = [
        ("pigz -dc", _pigz_decompress_command(in_path)),
        ("sort", _sort_command(ncores, sort_mem_pct)),
        ("uniq", ["uniq"]),
        ("bgzip", _bgzip_command(ncores)),
    ]

    logger.info(f"Starting deduplication pipeline: {ncores} cores, {sort_mem_pct}% memory per core")
    _pipeline_run(stages, out_path, env=_deterministic_env())
    logger.info(f"bgzip compression completed successfully for {out_path}")
    return str(out_path)
