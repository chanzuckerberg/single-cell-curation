## scATAC-seq asset Dataset Criteria

See [fragments file schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/drafts/fragments_file.md) for criteria a Dataset MUST meet in order to be eligible for scATAC-seq assets.

If a Dataset has a fragments file asset, it MAY have genome track assets. Otherwise, it MUST NOT have genome track assets.

## `uns` (Dataset Metadata)

<table><tbody>
    <tr>
      <th>Key</th>
      <td>peak_grouping</td>
    </tr>
    <tr>
      <th>Annotation</th>
      <td>Curator MAY annotate if the Dataset has a fragments file asset; otherwise, this key MUST NOT be present.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td>
          <code>str</code>. The value MUST match a key in <code>obs</code>. If annotated, genome track assets MUST be submitted. The following columns MUST NOT be specified:
      <ul>
        <li>assay_ontology_term_id</li>
        <li>cell_type_ontology_term_id</li>
        <li>development_stage_ontology_term_id</li>
        <li>disease_ontology_term_id</li>
        <li>organism_ontology_term_id</li>
        <li>self_reported_ethnicity_ontology_term_id</li>
        <li>sex_ontology_term_id</li>
        <li>tissue_ontology_term_id</li>
      </ul>
      Instead specify the corresponding Discover column such as <code>cell_type</code>.<br><br>
        </td>
    </tr>
</tbody></table>

## scATAC-seq Asset: Genome Track

If <code>uns['peak_grouping']</code> is annotated, there MUST be exactly one genome track asset submitted for each unique value in the obs column specified as determined by <code>anndata.obs.{peak_grouping_column}.unique()</code>. Otherwise, this MUST NOT be submitted.

Asset file specifications TBD based on the visualization solution. Accepting <a href="https://genome.ucsc.edu/goldenpath/help/bigWig.html">.bigWig format</a> is a requirement.
