from typing import List, Tuple


def merge_bed_ranges(ranges: List[Tuple]) -> List[Tuple]:
    """
    Merges bed-like ranges and returns non-overlapping ranges

    code adapted from http://www.genemine.org/gtftools.php

    :param ranges List[Tuple]: List of ranges, each range is a tuple of the form [chromosome, start, end, strand]

    :rtype List[Tuple]
    :return A list of non-overlapping ranges
    """

    ranges.sort(key=lambda x: (x[0], x[1]))
    merged = []
    n_range = len(ranges)
    if n_range == 1:
        merged = ranges
    elif n_range == 2:
        imerge = _neighbor_merge_ranges(ranges[0], ranges[1])
        for each in imerge:
            merged.append(each)
    else:
        i = 2
        imerge = _neighbor_merge_ranges(ranges[0], ranges[1])
        n_imerge = len(imerge)
        while n_imerge > 0:
            if n_imerge == 2:
                merged.append(imerge[0])

            imerge = _neighbor_merge_ranges(imerge[n_imerge - 1], ranges[i])
            n_imerge = len(imerge)
            if i == n_range - 1:
                for each in imerge:
                    merged.append(each)
                n_imerge = -1
            i += 1

    return merged


def _neighbor_merge_ranges(range1: Tuple, range2: Tuple) -> List[Tuple]:
    """
    Merges two neighboring ranges

    :param range1 Tuple: a tuple of the form [chromosome, start, end, strand]
    :param range2 Tuple: a tuple of the form [chromosome, start, end, strand]

    :rtype List[Tuple]
    :return A list of non-overlapping ranges
    """

    if range2[1] <= range1[2]:
        merged = [(range1[0], range1[1], max(range1[2], range2[2]), range1[3])]
    else:
        merged = [range1, range2]
    return merged


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
