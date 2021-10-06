import os
from typing import Union

from cellxgene_schema.validate import Validator
from cellxgene_schema.write_labels import AnnDataLabelAppender


def validate(h5ad_path: Union[str, bytes, os.PathLike], add_labels_file: str = None):

    """
    Entry point for validation.

    :param Union[str, bytes, os.PathLike] h5ad_path: Path to h5ad file to validate
    :param str add_labels_file: Path to new h5ad file with ontology/gene labels added

    :return True if successful validation, False otherwise
    :rtype bool
    """

    # Perform validation
    validator = Validator()
    validator.validate_adata(h5ad_path)

    # Stop if validation was unsuccessful
    if not validator.is_valid:
        return False

    if add_labels_file:
        writer = AnnDataLabelAppender(validator)
        writer.write_labels(add_labels_file)

        return validator.is_valid & writer.was_writing_successful

    return True
