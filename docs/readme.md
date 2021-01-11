## BICCN

### `BICCN_ontology_term_mappings.txt` 

Our efforts to integrate a newly born [brain ontology](https://scicrunch.org/scicrunch/interlex/view/ilx_0770272#relationships) adopted by BICCN have led to the creation and tracking of label mappings between the original dataset labels and the ontology terms. 

This table contains those mappings, with each row being a cell type mapping to an ontology term and the following columns: 

1. `original_columns`  what columns from `obs` in the original data were used for the mappings, if it's more than one column use `@` to separate the column names.
2. `original_values`  what values from those columns were mapped to the ontology term, if it's more than one column use `@` to separate the values.
3. `BICCN_ontology_term_id`  the mapped ontology id, if non copy the final value used.
4. `ontology_label`  the label assigned to that ontology.
5. `project_name` the name of dataset (using the first author, should be sufficient).
6. `doi` doi to the paper/preprint.

Ultimately our current goal to include this information on a `*.h5ad` file is to use it to add the following columns to the `obs` (cell metadata):

- `BICCN_class_label` (only if available) -- the most general level cell type definition from the original paper (so far GABAergic, Glutamergic, Non-Neuronal).
- `BICCN_subclass_label` -- intermediate level of cell type definition, these usually are the the labels they used throughout the paper and the ones mapping to the ontology (e.g. Endothelial, astrocyte, L5 ET, L5 IT).
- `BICCN_cluster_label`  (only if available)  -- the most granular cell type definition (e.g. L5 IT Tcap1, L5 IT Tcap2).
- `BICCN_ontology_term_id` -- the ontology term from from the table described above.
