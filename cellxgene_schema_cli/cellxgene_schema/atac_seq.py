import itertools
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import anndata as ad
import dask
import dask.dataframe as ddf
import dask.distributed as dd  # TODO: see if distributed mode can be avoided.
import h5py
import pandas as pd
import pysam
from dask import delayed
from dask.delayed import Delayed

from .ontology_parser import ONTOLOGY_PARSER
from .utils import is_ontological_descendant_of

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
allowed_chromosomes = list(set(itertools.chain(human_chromosome_by_length.keys(), mouse_chromosome_by_length.keys())))
allowed_chromosomes.sort()
chromosome_categories = pd.CategoricalDtype(categories=allowed_chromosomes)
column_types = {
    "chromosome": chromosome_categories,
    "start coordinate": int,
    "stop coordinate": int,
    "barcode": str,
    "read support": int,
}


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
    dask_cluster_config: Optional[dict] = None,
    override_write_algorithm: Optional[str] = None,
    output_file: Optional[str] = None,
) -> list[str]:
    """
    Validate the fragment against the anndata file and generate the index if the fragment is valid.

    :param str fragment_file: The fragment file to process
    :param str anndata_file: The anndata file to validate against
    :param bool generate_index: Whether to generate the index for the fragment
    :param dask_cluster_config: dask cluster configuration parameters passed to dask.distributed.LocalCluster
    :param override_write_algorithm: Override the write algorithm used to write the bgzip file. Options are "pysam"
    and "cli"
    :param output_file: The output file to write the bgzip file to. If not provided, the output file will be the same

    """
    with tempfile.TemporaryDirectory() as tempdir:

        # configure the dask cluster
        _dask_cluster_config = dict(dashboard_address=None)
        logging.getLogger("distributed").setLevel(logging.ERROR)
        if dask_cluster_config:
            _dask_cluster_config.update(dask_cluster_config)

        # start the dask cluster and client
        with dd.LocalCluster(**_dask_cluster_config) as cluster, dd.Client(cluster):
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

            # slow checks
            errors = validate_anndata_with_fragment(parquet_file, anndata_file)
            if errors:
                return errors
            else:
                logger.info("Fragment and Anndata file are valid")

            # generate the index
            if generate_index:
                logger.info(f"Sorting fragment and generating index for {fragment_file}")
                index_fragment(fragment_file, parquet_file, tempdir, override_write_algorithm, output_file)
        logger.debug("cleaning up")
    return []


def convert_to_parquet(fragment_file: str, tempdir: str) -> str:
    # unzip the fragment. Subprocess is used here because gzip on the cli uses less memory with comparable
    # speed to the python gzip library.
    unzipped_file = Path(tempdir) / Path(fragment_file).name.rsplit(".", 1)[0]
    logger.info(f"Unzipping {fragment_file}")
    with open(unzipped_file, "wb") as fp:
        subprocess.run(["gunzip", "-c", fragment_file], stdout=fp, check=True)

    # convert the fragment to a parquet file
    logger.info(f"Converting {fragment_file} to parquet")
    parquet_file_path = Path(tempdir) / Path(fragment_file).name.replace(".gz", ".parquet")
    try:
        ddf.read_csv(
            unzipped_file, sep="\t", names=column_ordering, dtype=column_types, keep_default_na=False
        ).to_parquet(parquet_file_path, partition_on=["chromosome"], compute=True)
    finally:
        # remove the unzipped file
        logger.debug(f"Removing {unzipped_file}")
        unzipped_file.unlink()
    parquet_file = str(parquet_file_path)
    return parquet_file


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


def validate_anndata_with_fragment(parquet_file: str, anndata_file: str) -> list[str]:
    errors = [
        validate_fragment_start_coordinate_greater_than_0(parquet_file),
        validate_fragment_barcode_in_adata_index(parquet_file, anndata_file),
        validate_fragment_stop_coordinate_within_chromosome(parquet_file, anndata_file),
        validate_fragment_stop_greater_than_start_coordinate(parquet_file),
        validate_fragment_read_support(parquet_file),
        validate_fragment_no_duplicate_rows(parquet_file),
    ]
    return report_errors("Errors found in Fragment and/or Anndata file", errors)


def validate_fragment_no_duplicate_rows(parquet_file: str) -> Optional[str]:
    df = ddf.read_parquet(parquet_file)
    if len(df.drop_duplicates()) != len(df):
        return "Fragment file has duplicate rows."


def validate_fragment_start_coordinate_greater_than_0(parquet_file: str) -> Optional[str]:
    df = ddf.read_parquet(parquet_file, columns=["start coordinate"])
    series = df["start coordinate"] > 0
    if not series.all().compute():
        return "Start coordinate must be greater than 0."


def validate_fragment_barcode_in_adata_index(parquet_file: str, anndata_file: str) -> Optional[str]:
    df = ddf.read_parquet(parquet_file, columns=["barcode"])
    with h5py.File(anndata_file) as f:
        obs = ad.io.read_elem(f["obs"])
    barcode = set(df.groupby(by="barcode").count().compute().index)
    if set(obs.index) != barcode:
        return "Barcodes don't match anndata.obs.index"


def validate_fragment_stop_greater_than_start_coordinate(parquet_file: str) -> Optional[str]:
    df = ddf.read_parquet(parquet_file, columns=["start coordinate", "stop coordinate"])
    series = df["stop coordinate"] > df["start coordinate"]
    if not series.all().compute():
        return "Stop coordinate must be greater than start coordinate."


def validate_fragment_stop_coordinate_within_chromosome(parquet_file: str, anndata_file: str) -> Optional[str]:
    # check that the stop coordinate is within the length of the chromosome
    with h5py.File(anndata_file) as f:
        organism_ontology_term_id = ad.io.read_elem(f["obs"])["organism_ontology_term_id"].unique().astype(str)[0]
    df = ddf.read_parquet(parquet_file, columns=["chromosome", "stop coordinate"])

    # the chromosomes in the fragment must match the chromosomes for that organism
    chromosome_length_table = organism_ontology_term_id_by_chromosome_length_table[organism_ontology_term_id]
    mismatched_chromosomes = set(df["chromosome"].unique().compute()) - chromosome_length_table.keys()
    if mismatched_chromosomes:
        return f"Chromosomes in the fragment do not match the organism({organism_ontology_term_id}).\n" + "\t\n".join(
            mismatched_chromosomes
        )
    df["chromosome_length"] = df["chromosome"].map(chromosome_length_table, meta=int).astype(int)
    df = df["stop coordinate"] <= df["chromosome_length"]
    if not df.all().compute():
        return "Stop coordinate must be less than the chromosome length."


def validate_fragment_read_support(parquet_file: str) -> Optional[str]:
    # check that the read support is greater than 0
    df = ddf.read_parquet(
        parquet_file, columns=["read support"], filters=[("read support", "<=", 0)], keep_default_na=False
    )
    if len(df.compute()) != 0:
        return "Read support must be greater than 0."


def validate_anndata_is_primary_data(anndata_file: str) -> Optional[str]:
    with h5py.File(anndata_file) as f:
        is_primary_data = ad.io.read_elem(f["obs"])["is_primary_data"]
    if not is_primary_data.all():
        return "Anndata.obs.is_primary_data must all be True."


def validate_anndata_organism_ontology_term_id(anndata_file: str) -> Optional[str]:
    with h5py.File(anndata_file) as f:
        organism_ontology_term_ids = ad.io.read_elem(f["obs"])["organism_ontology_term_id"].unique().astype(str)
    if organism_ontology_term_ids.size > 1:
        error_message = (
            "Anndata.obs.organism_ontology_term_id must have exactly 1 unique value. Found the following values:\n"
        ) + "\n\t".join(organism_ontology_term_ids)
        return error_message
    organism_ontology_term_id = organism_ontology_term_ids[0]
    allowed_terms = [*organism_ontology_term_id_by_chromosome_length_table.keys()]
    if organism_ontology_term_id not in allowed_terms:
        return f"Anndata.obs.organism_ontology_term_id must be one of {allowed_terms}. Got {organism_ontology_term_id}."


def detect_chromosomes(parquet_file: str) -> list[str]:
    logger.info("detecting chromosomes")
    df = ddf.read_parquet(parquet_file, columns=["chromosome"]).drop_duplicates()
    df = df["chromosome"].compute()
    return df.tolist()


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
    parquet_file: str,
    tempdir: str,
    override_write_algorithm: Optional[str] = None,
    output_file: Optional[str] = None,
):
    # sort the fragment by chromosome, start coordinate, and stop coordinate, then compress it with bgzip
    bgzip_output_file = get_output_file(fragment_file, output_file)
    bgzip_output_path = Path(bgzip_output_file)
    bgzip_output_path.unlink(missing_ok=True)
    bgzip_output_path.touch()
    bgzip_write_lock = dd.Lock()  # lock to preserver write order

    if override_write_algorithm:
        write_algorithm = write_algorithm_by_callable[override_write_algorithm]
    elif not shutil.which("bgzip"):  # check if bgzip cli is installed
        logger.warning("bgzip is not installed, using slower pysam implementation")
        write_algorithm = write_algorithm_by_callable["pysam"]
    else:
        write_algorithm = write_algorithm_by_callable["cli"]

    chromosomes = detect_chromosomes(parquet_file)
    jobs = prepare_fragment(chromosomes, parquet_file, bgzip_output_file, tempdir, bgzip_write_lock, write_algorithm)
    # limit calls to dask.compute to improve performance. The number of jobs to run at once is determined by the
    # step variable. If we run all the jobs in the same call to dask.compute, the local cluster hangs.
    # TODO: investigate why
    step = 4
    # print the progress of the jobs
    for i in range(0, len(jobs), step):
        dask.compute(jobs[i : i + step])

    logger.info(f"Fragment sorted and compressed: {bgzip_output_file}")
    #
    pysam.tabix_index(bgzip_output_file, preset="bed", force=True)
    tabix_output_file = bgzip_output_file + ".tbi"
    logger.info(f"Index file generated: {tabix_output_file}")


@delayed
def sort_fragment(parquet_file: str, write_path: str, chromosome: str) -> Path:
    temp_data = Path(write_path) / f"temp_{chromosome}.tsv.gz"
    df = ddf.read_parquet(parquet_file, filters=[("chromosome", "==", chromosome)])
    df = df[column_ordering]
    df = df.sort_values(["start coordinate", "stop coordinate"], ascending=True)

    df.to_csv(temp_data, sep="\t", index=False, header=False, mode="w", single_file=True)
    return temp_data


@delayed
def write_bgzip_pysam(input_file: str, bgzip_output_file: str, write_lock: dd.Lock):
    with write_lock, pysam.libcbgzf.BGZFile(bgzip_output_file, mode="ab") as f_out, open(input_file, "rb") as f_in:
        while data := f_in.read(2**20):
            f_out.write(data)


@delayed
def write_bgzip_cli(input_file: str, bgzip_output_file: str, write_lock: dd.Lock):
    with write_lock, open(input_file, "rb") as fin, open(bgzip_output_file, "ab") as fout:
        subprocess.run(["bgzip", "--threads=8", "-c"], stdin=fin, stdout=fout, check=True)


write_algorithm_by_callable = {"pysam": write_bgzip_pysam, "cli": write_bgzip_cli}


def prepare_fragment(
    chromosomes: list[str],
    parquet_file: str,
    bgzip_output_file: str,
    tempdir: str,
    write_lock: dd.Lock,
    write_algorithm: callable,
) -> list[Delayed]:
    """
    The sorting and writing of the fragment is done in parallel for each chromosome. Because of this the write order of
    the chromosomes may not be preserved. The chromosomes will all be stored in contiguous blocks in the bgzip file, and
    and sorted by start and stop coordinate within each chromosome. This slightly different form what pysam.tabix
    expects which is sorted by chromosome, start coordinate, and stop coordinate, but it is still compatible.

    :param chromosomes:
    :param parquet_file:
    :param bgzip_output_file:
    :param tempdir:
    :param write_lock:
    :param write_algorithm:
    :return:
    """
    jobs = []
    for chromosome in chromosomes:
        temp_data = sort_fragment(parquet_file, tempdir, chromosome)
        jobs.append(write_algorithm(temp_data, bgzip_output_file, write_lock))
    return jobs
