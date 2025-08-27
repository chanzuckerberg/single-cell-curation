import gzip
import logging
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Optional

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


def line_counts_match(file1: str, file2: str) -> None:
    """
    Validate that the line counts of two compressed files match.
    """
    logger.info("Validating line count between original and output files")
    with ThreadPoolExecutor(max_workers=2) as executor:
        original_future = executor.submit(count_lines_in_compressed_file, file1)
        output_future = executor.submit(count_lines_in_compressed_file, file2)

        original_line_count = original_future.result()
        output_line_count = output_future.result()

    if original_line_count != output_line_count:
        error_msg = f"Line count validation failed: original file has {original_line_count} lines, output file has {output_line_count} lines"
        logger.error(error_msg)
        raise ValueError(error_msg)

    logger.info(f"Line count validation passed: {original_line_count} lines in both original and output files")


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
    fragment_is_prepared: bool = False,
) -> list[str]:
    """
    Validate the fragment against the anndata file and generate the index if the fragment is valid.

    :param str fragment_file: The fragment file to process
    :param str anndata_file: The anndata file to validate against
    :param bool generate_index: Whether to generate the index for the fragment
    :param override_write_algorithm: Override the write algorithm used to write the bgzip file. Options are "pysam"
    and "cli"
    :param output_file: The output file to write the bgzip file to. If not provided, the output file will be the same
    :param fragment_is_prepared: If True, skip the preparation step and assume the fragment is already sorted and bgzipped.
    :return: A list of error messages. If the list is empty, the fragment is valid.
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
        if not fragment_is_prepared:
            logger.info("Preparing fragment for indexing")
            prepared_fragment_file = prepare_fragment(fragment_file, output_file)
        else:
            logger.info("Fragment is already prepared, skipping preparation step")
            prepared_fragment_file = fragment_file

        logger.info(f"Sorting fragment and generating index for {prepared_fragment_file}")
        index_fragment(prepared_fragment_file)

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


def index_fragment(bgzip_fragment_file: str) -> None:
    pysam.tabix_index(bgzip_fragment_file, preset="bed", force=True)
    tabix_output_file = bgzip_fragment_file + ".tbi"
    logger.info(f"Index file generated: {tabix_output_file}")


SORT_MEMORY_PERCENTAGE = (
    80  # percentage of memory to use for sort command. This is the default used by linux sort command.
)


def prepare_fragment(
    fragment_file: str,
    bgzip_output_file: str,
) -> str:
    """
    The sorting and writing of the fragment is done for each chromosome. Because of this the write order of
    the chromosomes may not be preserved. The chromosomes will all be stored in contiguous blocks in the bgzip file, and
    and sorted by start and stop coordinate within each chromosome.
    :param fragment_file: The fragment file to process. This is a gzipped compressed fragment file.
    :param bgzip_output_file: The output file to write the bgzip file to.
    :return: The path to the bgzip file.
    """

    check_external_requirements()

    # remove existing output file if it exists
    bgzip_output_file = get_output_file(fragment_file, bgzip_output_file)
    bgzip_output_path = Path(bgzip_output_file)
    bgzip_output_path.unlink(missing_ok=True)

    logger.info(f"Fragment sorted and compressed: {bgzip_output_file}")
    num_cores = os.cpu_count()
    with open(bgzip_output_file, "wb") as out_f:
        gzip_proc = subprocess.Popen(get_gzip_command(fragment_file), stdout=subprocess.PIPE)
        sort_proc = subprocess.Popen(
            get_sort_command(num_cores),
            stdin=gzip_proc.stdout,
            stdout=subprocess.PIPE,
            env={**os.environ, "LC_ALL": "C"},
        )
        gzip_proc.stdout.close()
        bgzip_proc = subprocess.Popen(get_bgzip_command(num_cores), stdin=sort_proc.stdout, stdout=out_f)
        sort_proc.stdout.close()
        bgzip_proc.wait()
    if bgzip_proc.returncode != 0:
        raise RuntimeError(f"bgzip compression failed with error code {bgzip_proc.returncode}")
    line_counts_match(fragment_file, bgzip_output_file)
    logger.info(f"bgzip compression completed successfully for {bgzip_output_file}")
    return bgzip_output_file


def get_sort_command(num_cores: int, sort_memory_percent: int = SORT_MEMORY_PERCENTAGE) -> list[str]:
    """
    Get the sort command with the appropriate parameters for parallelization and memory usage.

    :param num_cores: Number of CPU cores available.
    :param sort_memory_percent: Total percentage of system memory to allocate for sorting.
    :return: List of command arguments for the sort command.
    """
    sort_memory = calculate_sort_memory(num_cores, sort_memory_percent)
    return [
        "sort",
        "-t",
        "\t",
        "-k1,1",
        "-k2,2n",
        "-k3,3n",
        "-k4,4",
        f"-S{sort_memory}%",
        f"--parallel={num_cores}",
        f'--compress-program="pigz -p {num_cores}"',
    ]


def get_bgzip_command(num_cores: int) -> list[str]:
    """
    Get the bgzip command with the appropriate parameters for parallelization.

    :param num_cores: Number of CPU cores available.
    :return: List of command arguments for the bgzip command.
    """
    return ["bgzip", f"--threads={num_cores}", "-c"]


def get_gzip_command(file_name: str) -> list[str]:
    """
    Get the gzip command.

    :param file_name: The file to compress.
    :return: List of command arguments for the gzip command.
    """
    return ["gzip", "-dc", file_name]


def calculate_sort_memory(num_cores: int, sort_memory_percent: int) -> int:
    """
    Calculate the memory percentage to allocate per sort thread.

    :param num_cores: Number of CPU cores available.
    :param sort_memory_percent: Total percentage of system memory to allocate for sorting.
    :return: Percentage of memory to allocate per sort thread.
    """
    return max(sort_memory_percent // num_cores, 1)  # ensure at least 1% memory per core


def check_external_requirements() -> None:
    if shutil.which("sort") is None:
        raise RuntimeError(
            "The 'sort' command is not installed or not found in PATH. It is required to sort the fragment file."
        )
    if shutil.which("bgzip") is None:
        raise RuntimeError(
            "The 'bgzip' command is not installed or not found in PATH. It is required to compress the fragment file."
        )
    if shutil.which("pigz") is None:
        raise RuntimeError(
            "The 'pigz' command is not installed or not found in PATH. It is required to compress the fragment file."
        )
    if shutil.which("uniq") is None:
        raise RuntimeError(
            "The 'uniq' command is not installed or not found in PATH. It is required to compress the fragment file."
        )


def deduplicate_fragment_rows(
    fragment_file_name: str, output_file_name: str = None, sort_memory_percent: int = SORT_MEMORY_PERCENTAGE
) -> str:
    """
    Deduplicate rows in a fragment file by sorting and using the uniq command, then compress the result with bgzip.

    This function decompresses the input fragment file, sorts it by chromosome and coordinates using GNU sort
    (with parallelization and pigz compression), removes duplicate rows, and writes the output as a bgzipped file.
    The amount of memory allocated to sort is controlled by sort_memory_percent, divided among CPU cores.

    :param fragment_file_name: Path to the input gzipped fragment file.
    :param output_file_name: Path for the output deduplicated bgzipped file. If None, a default name is generated.
    :param sort_memory_percent: Percentage of system memory to allocate per sort thread (default: 80).
    :return: Path to the deduplicated bgzipped output file.
    :raises RuntimeError: If the bgzip compression process fails.
    """

    check_external_requirements()
    logger.info(f"deduplicating fragment file: {fragment_file_name}")

    output_file_name = (
        Path(output_file_name)
        if output_file_name
        else (Path(fragment_file_name).stem.replace(".tsv", "") + "_dedup.tsv.bgz")
    )

    num_cores = os.cpu_count()
    # Build the pipeline: gzip -dc | sort | uniq | bgzip
    gzip_cmd = get_gzip_command(fragment_file_name)
    sort_cmd = get_sort_command(num_cores, sort_memory_percent)
    uniq_cmd = ["uniq"]
    bgzip_cmd = get_bgzip_command(num_cores)

    with open(output_file_name, "wb") as outfile:
        # Start gzip process
        gzip_proc = subprocess.Popen(gzip_cmd, stdout=subprocess.PIPE)
        # Start sort process
        sort_proc = subprocess.Popen(
            sort_cmd,
            stdin=gzip_proc.stdout,
            stdout=subprocess.PIPE,
            env={**os.environ, "LC_ALL": "C"},
        )
        gzip_proc.stdout.close()  # Allow gzip_proc to receive a SIGPIPE if sort_proc exits.
        # Start uniq process
        uniq_proc = subprocess.Popen(uniq_cmd, stdin=sort_proc.stdout, stdout=subprocess.PIPE)
        sort_proc.stdout.close()
        # Start bgzip process
        bgzip_proc = subprocess.Popen(bgzip_cmd, stdin=uniq_proc.stdout, stdout=outfile)
        uniq_proc.stdout.close()
        bgzip_proc.wait()
        if bgzip_proc.returncode != 0:
            raise RuntimeError(f"bgzip compression failed with error code {bgzip_proc.returncode}")
        # Wait for the rest of the pipeline to finish and check for errors
        uniq_proc.wait()
        sort_proc.wait()
        gzip_proc.wait()
    logger.info(f"bgzip compression completed successfully for {output_file_name}")
    return str(output_file_name)
