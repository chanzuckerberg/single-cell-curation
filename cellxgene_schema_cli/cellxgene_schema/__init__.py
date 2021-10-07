import logging
import os
from datetime import datetime
from typing import Union

from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender

logger = logging.getLogger(__name__)


def validate(h5ad_path: Union[str, bytes, os.PathLike], add_labels_file: str = None, verbose: bool = False) -> (
        bool, list):
    """
    Entry point for validation.

    :param Union[str, bytes, os.PathLike] h5ad_path: Path to h5ad file to validate
    :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

    :return (True, []) if successful validation, (False, [list_of_errors] otherwise
    :rtype tuple
    """

    # Perform validation
    start = datetime.now()
    if verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    validator = Validator()
    validator.validate_adata(h5ad_path)
    logger.debug(f"Validation complete in {datetime.now() - start} seconds with status is_valid={validator.is_valid}")

    # Stop if validation was unsuccessful
    if not validator.is_valid:
        return False, validator.errors

    if add_labels_file:
        label_start = datetime.now()
        writer = AnnDataLabelAppender(validator)
        writer.write_labels(add_labels_file)
        logger.debug(f"H5AD label writing complete in {datetime.now() - label_start}, was_writing_successful: {writer.was_writing_successful}") # noqa E501

        return validator.is_valid & writer.was_writing_successful, validator.errors

    return True
