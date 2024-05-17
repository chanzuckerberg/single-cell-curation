from typing import List

import anndata as ad
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from tests.fixtures.examples_validate import h5ad_valid


class ConfiguredBaseModel(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)


class UnsModel(ConfiguredBaseModel):
    title: str
    batch_condition: List[str]
    # colors: ndarray


class AnndataModel(ConfiguredBaseModel):
    adata: ad.AnnData

    @model_validator(mode="before")
    @classmethod
    def _validate_encoding_version(cls, adata: ad.AnnData):
        import h5py

        with h5py.File(adata.filename, "r") as f:
            encoding_dict = dict(f.attrs)
            encoding_version = encoding_dict.get("encoding-version")
            if encoding_version != "0.1.0":
                raise ValueError("The h5ad artifact was generated with an AnnData version different from 0.8.0.")
        return adata


class CXGModel(ConfiguredBaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True, from_attributes=True)
    adata: AnndataModel
    uns: UnsModel
    # X: Union[np.ndarray, sparse.spmatrix, ArrayView, SparseDataset, None]
    # obs: Union[pd.DataFrame, Mapping[str, Iterable[Any]], None]
    # var: Union[pd.DataFrame, Mapping[str, Iterable[Any]], None]
    #
    # obsm: Union[np.ndarray, Mapping[str, Union[Sequence[Any], np.ndarray]], None]
    # varm: Union[np.ndarray, Mapping[str, Union[Sequence[Any], np.ndarray]], None]
    # layers: Union[Layers , LayersView]
    # raw: Union[Mapping[str, Any], None, Raw]
    # filename: Union[PathLike, None]
    # obsp: Union[np.ndarray, Mapping[str, Union[Sequence[Any], np.ndarray]], None]
    # varp: Union[np.ndarray, Mapping[str, Union[Sequence[Any], np.ndarray]], None]


data = ad.read_h5ad(h5ad_valid, backed="r")
v = AnndataModel.model_validate(data)


# we have a lot of custom types like overloadedDict
# unable to validate the memeber of the anndata class


class EmbeddingSchema(BaseModel):
    title: str
    batch_condition: List[str]
    default_embedding: str
    X_approximate_distribution: str

    @model_validator(mode="after")
    def validate_forbidden_keys(cls, values):
        forbidden_keys = [
            "schema_version",
            "citation",
            "schema_reference",
            "X_normalization",
            "default_field",
            "layer_descriptions",
            "tags",
            "version",
            "contributors",
            "preprint_doi",
            "project_description",
            "project_links",
            "project_name",
            "publication_doi",
        ]
        for key in forbidden_keys:
            if key in values:
                raise ValueError(f"Forbidden key found: {key}")
        return values

    @field_validator("X_approximate_distribution")
    def validate_X_approximate_distribution(cls, value):
        if value not in ["count", "normal"]:
            raise ValueError("X_approximate_distribution must be 'count' or 'normal'")
        return value


class ColorsSchema(BaseModel):
    category: str = Field(pattern=r"^.*(_colors)|(_ontology_term_id_colors)$")
    colors: List[str]

    @model_validator(mode="after")
    def validate_forbidden_keys(cls, values):
        forbidden_keys = [
            "assay_colors",
            "cell_type_colors",
            "development_stage_colors",
            "disease_colors",
            "organism_colors",
            "self_reported_ethnicity_colors",
            "sex_colors",
            "tissue_colors",
        ]
        for key in forbidden_keys:
            if key in values:
                raise ValueError(f"Forbidden key found: {key}")
        return values

    @model_validator(mode="after")
    def validate_colors(cls, values):
        # Custom validation logic for colors
        # Implement the "obs_category_mapping" check here
        return values


# Test the models
if __name__ == "__main__":
    # Example input data
    example_data_1 = {
        "title": "Example Title",
        "batch_condition": ["condition1", "condition2"],
        "default_embedding": "embedding_key",
        "X_approximate_distribution": "count",
        "schema_version": "1.0",  # Forbidden key
    }
    example_data_2 = {
        "category": "assay_colors",
        "colors": ["#ffffff", "#000000"],
        "assay_colors": ["#ffffff", "#000000"],  # Forbidden key
    }

    # Create and validate EmbeddingSchema instance
    embedding_instance = EmbeddingSchema(**example_data_1)
    print(embedding_instance.dict())

    # Create and validate ColorsSchema instance
    colors_instance = ColorsSchema(**example_data_2)
    print(colors_instance.dict())
