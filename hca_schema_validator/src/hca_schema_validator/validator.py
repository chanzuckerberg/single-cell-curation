"""HCA Validator - extends cellxgene Validator with HCA-specific rules."""

from cellxgene_schema.validate import Validator


class HCAValidator(Validator):
    """
    HCA-specific validator extending cellxgene schema validation.
    
    Currently a passthrough to the base Validator.
    Future enhancements will add HCA-specific validation rules.
    """
    
    def __init__(self, ignore_labels=False):
        """
        Initialize HCA validator.
        
        Args:
            ignore_labels: If True, skip label validation
        """
        super().__init__(ignore_labels=ignore_labels)
