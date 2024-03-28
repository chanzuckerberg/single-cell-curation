import json
import logging
from urllib.request import urlopen

from cellxgene_ontology_guide.ontology_parser import OntologyParser
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)
    ontology_parser = OntologyParser()
    datasets = json.loads(urlopen("https://api.cellxgene.cziscience.com/dp/v1/datasets/index").read().decode("utf-8"))
    with logging_redirect_tqdm():
        for dataset in tqdm(datasets):
            logger_dataset = logging.getLogger(f"dataset_{dataset['id']}")
            list_of_terms_and_labels = [
                *dataset.get("cell_type", []),
                *dataset.get("development_stage", []),
                *dataset.get("disease", []),
                *dataset.get("organism", []),
                *dataset.get("self_reported_ethnicity", []),
                *dataset.get("sex", []),
                *dataset.get("tissue", []),
            ]
            for term_and_label in list_of_terms_and_labels:
                if "," in term_and_label["ontology_term_id"]:
                    terms = term_and_label["ontology_term_id"].split(",")
                    labels = [i.lstrip() for i in term_and_label["label"].split(",")]
                else:
                    terms = [term_and_label["ontology_term_id"]]
                    labels = [term_and_label["label"]]
                for term in terms:
                    cog_term_label = ontology_parser.get_term_label(term)
                    if cog_term_label not in labels:
                        logger_dataset.error(f"index_term_label='{labels}' != {cog_term_label=}")
                    else:
                        labels.remove(cog_term_label)
                if labels:
                    logger_dataset.error("extra labels found: " + ", ".join(labels))


if __name__ == "__main__":
    main()
