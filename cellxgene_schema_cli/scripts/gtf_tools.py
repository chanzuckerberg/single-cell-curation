def _get_features(gtf_line: str) -> dict:
    """
    Parses the features found in column 8 of GTF, returns a dict with keys as feature names and values as the feature
    values

    :param str gtf_line: a line from a GTF file

    :rtype dict
    :Return a dictionary with keys and values from column 8 of GTF; key-value example: gene_id-ENSG0000001
    """

    return_features = {}

    for feature in gtf_line[8].split(";"):
        if len(feature) > 1:
            feature = feature.strip().split(" ")
            feature_name = feature[0]
            feature_value = feature[1]
            return_features[feature_name] = feature_value.replace('"', "")

    return return_features
