## scATAC-seq assay types

<i>paired assay</i>: any descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010891"><code>"EFO:0010891"</code></a> for <i>scATAC-seq</i> that is also a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008913"><code>"EFO:0008913"</code></a> for <i>single-cell RNA sequencing</i>

<i>unpaired assay</i>: <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0010891"><code>"EFO:0010891"</code></a> for <i>scATAC-seq</i> or its descendants that is not a descendant of <a href="https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=EFO%3A0008913"><code>"EFO:0008913"</code></a> for <i>single-cell RNA sequencing</i>

## Fragment File Dataset Criteria

A Dataset MUST meet each of the following criteria in order to be eligible for an attached Fragment File:
* the <code>obs['assay_ontology_term_id']</code> values MUST all be <i>paired assays</i> or MUST all be <i>unpaired assays</i>
* the <code>obs['is_primary_data']</code> values MUST be all `True`
* the <code>var['feature_reference']</code> values MUST include one of <code>"NCBITaxon:9606"</code> for <i>Homo sapiens</i> or <code>"NCBITaxon:10090"</code> for <i>Mus musculus</i>, but not both. The value that is present will determine the appropriate Chromosome Table for standards.

If the <code>obs['assay_ontology_term_id']</code> values are all <i>paired assays</i> then a fragment file MAY be attached to the Dataset.

If the <code>obs['assay_ontology_term_id']</code> values are all <i>unpaired assays</i> then a fragment file MUST be attached to the Dataset.

## Fragment File

This MUST be a gzipped tab-separated values (TSV) file.

The curator MUST annotate the following header-less columns. Additional columns and header lines beginning with `#` MUST NOT be included.

### first column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. This MUST be the reference genome chromosome the fragment is located on. The value MUST be one of the values from the <code>Chromosome</code> column in the appropriate Chromosome Table.
        </td>
    </tr>
</tbody></table>
<br>

### second column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be the 0-based start coordinate of the fragment.
        </td>
    </tr>
</tbody></table>
<br>

### third column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be the 0-based end coordinate of the fragment. The end position is exclusive, so represents the position immediately following the fragment interval. The value MUST be greater than the start coordinate specified in the second column and less than or equal to the <code>Length</code> of the <code>Chromosome</code> specified in the first column, as specified in the appropriate Chromosome Table.
        </td>
    </tr>
</tbody></table>
<br>

### fourth column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>str</code>. This MUST be the cell identifier. The value MUST be found in the <code>obs</code> index of the associated Dataset. Every <code>obs</code> index value of the associated Dataset MUST appear at least once in this column of the fragment file.
        </td>
    </tr>
</tbody></table>
<br>

### fifth column

<table><tbody>
    <tr>
      <th>Annotator</th>
      <td>Curator MUST annotate.</td>
    </tr>
    <tr>
      <th>Value</th>
        <td><code>int</code>. This MUST be the total number of read pairs associated with this fragment. The value MUST be <code>1</code> or greater.
        </td>
    </tr>
</tbody></table>
<br>

## Fragment File index

CELLxGENE Discover MUST generate a <a href="https://www.htslib.org/doc/tabix.html">tabix</a> index of the fragment intervals from the fragment file. The file name MUST be the name of the corresponding fragment file appended with `.tbi`.

## Chromosome Tables

As determined by the reference assembly used by the gene annotation versions pinned for this version of the schema. Only chromosomes or scaffolds that have at least one gene feature present are included.

### <a href="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz">human (GRCh38.p14)</a>

<table>
  <thead>
  <tr>
  <th>Chromosome</th>
  <th>Length</th>
  </tr>
  </thead>
  <tbody>
    <tr>
        <td>chr1</td>
        <td>248956422</td>
    </tr>
    <tr>
        <td>chr2</td>
        <td>242193529</td>
    </tr>
    <tr>
        <td>chr3</td>
        <td>198295559</td>
    </tr>
    <tr>
        <td>chr4</td>
        <td>190214555</td>
    </tr>
    <tr>
        <td>chr5</td>
        <td>181538259</td>
    </tr>
    <tr>
        <td>chr6</td>
        <td>170805979</td>
    </tr>
    <tr>
        <td>chr7</td>
        <td>159345973</td>
    </tr>
    <tr>
        <td>chr8</td>
        <td>145138636</td>
    </tr>
    <tr>
        <td>chr9</td>
        <td>138394717</td>
    </tr>
    <tr>
        <td>chr10</td>
        <td>133797422</td>
    </tr>
    <tr>
        <td>chr11</td>
        <td>135086622</td>
    </tr>
    <tr>
        <td>chr12</td>
        <td>133275309</td>
    </tr>
    <tr>
        <td>chr13</td>
        <td>114364328</td>
    </tr>
    <tr>
        <td>chr14</td>
        <td>107043718</td>
    </tr>
    <tr>
        <td>chr15</td>
        <td>101991189</td>
    </tr>
    <tr>
        <td>chr16</td>
        <td>90338345</td>
    </tr>
    <tr>
        <td>chr17</td>
        <td>83257441</td>
    </tr>
    <tr>
        <td>chr18</td>
        <td>80373285</td>
    </tr>
    <tr>
        <td>chr19</td>
        <td>58617616</td>
    </tr>
    <tr>
        <td>chr20</td>
        <td>64444167</td>
    </tr>
    <tr>
        <td>chr21</td>
        <td>46709983</td>
    </tr>
    <tr>
        <td>chr22</td>
        <td>50818468</td>
    </tr>
    <tr>
        <td>chrX</td>
        <td>156040895</td>
    </tr>
    <tr>
        <td>chrY</td>
        <td>57227415</td>
    </tr>
    <tr>
        <td>chrM</td>
        <td>16569</td>
    </tr>
    <tr>
        <td>GL000009.2</td>
        <td>201709</td>
    </tr>
    <tr>
        <td>GL000194.1</td>
        <td>191469</td>
    </tr>
    <tr>
        <td>GL000195.1</td>
        <td>182896</td>
    </tr>
    <tr>
        <td>GL000205.2</td>
        <td>185591</td>
    </tr>
    <tr>
        <td>GL000213.1</td>
        <td>164239</td>
    </tr>
    <tr>
        <td>GL000216.2</td>
        <td>176608</td>
    </tr>
    <tr>
        <td>GL000218.1</td>
        <td>161147</td>
    </tr>
    <tr>
        <td>GL000219.1</td>
        <td>179198</td>
    </tr>
    <tr>
        <td>GL000220.1</td>
        <td>161802</td>
    </tr>
    <tr>
        <td>GL000225.1</td>
        <td>211173</td>
    </tr>
    <tr>
        <td>KI270442.1</td>
        <td>392061</td>
    </tr>
    <tr>
        <td>KI270711.1</td>
        <td>42210</td>
    </tr>
    <tr>
        <td>KI270713.1</td>
        <td>40745</td>
    </tr>
    <tr>
        <td>KI270721.1</td>
        <td>100316</td>
    </tr>
    <tr>
        <td>KI270726.1</td>
        <td>43739</td>
    </tr>
    <tr>
        <td>KI270727.1</td>
        <td>448248</td>
    </tr>
    <tr>
        <td>KI270728.1</td>
        <td>1872759</td>
    </tr>
    <tr>
        <td>KI270731.1</td>
        <td>150754</td>
    </tr>
    <tr>
        <td>KI270733.1</td>
        <td>179772</td>
    </tr>
    <tr>
        <td>KI270734.1</td>
        <td>165050</td>
    </tr>
    <tr>
        <td>KI270744.1</td>
        <td>168472</td>
    </tr>
    <tr>
        <td>KI270750.1</td>
        <td>148850</td>
    </tr>
</tbody></table>

### <a href="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/GRCm39.primary_assembly.genome.fa.gz">mouse (GRCm39)</a>

<table>
  <thead>
  <tr>
  <th>Chromosome</th>
  <th>Length</th>
  </tr>
  </thead>
  <tbody>
    <tr>
        <td>chr1</td>
        <td>195154279</td>
    </tr>
    <tr>
        <td>chr2</td>
        <td>181755017</td>
    </tr>
    <tr>
        <td>chr3</td>
        <td>159745316</td>
    </tr>
    <tr>
        <td>chr4</td>
        <td>156860686</td>
    </tr>
    <tr>
        <td>chr5</td>
        <td>151758149</td>
    </tr>
    <tr>
        <td>chr6</td>
        <td>149588044</td>
    </tr>
    <tr>
        <td>chr7</td>
        <td>144995196</td>
    </tr>
    <tr>
        <td>chr8</td>
        <td>130127694</td>
    </tr>
    <tr>
        <td>chr9</td>
        <td>124359700</td>
    </tr>
    <tr>
        <td>chr10</td>
        <td>130530862</td>
    </tr>
    <tr>
        <td>chr11</td>
        <td>121973369</td>
    </tr>
    <tr>
        <td>chr12</td>
        <td>120092757</td>
    </tr>
    <tr>
        <td>chr13</td>
        <td>120883175</td>
    </tr>
    <tr>
        <td>chr14</td>
        <td>125139656</td>
    </tr>
    <tr>
        <td>chr15</td>
        <td>104073951</td>
    </tr>
    <tr>
        <td>chr16</td>
        <td>98008968</td>
    </tr>
    <tr>
        <td>chr17</td>
        <td>95294699</td>
    </tr>
    <tr>
        <td>chr18</td>
        <td>90720763</td>
    </tr>
    <tr>
        <td>chr19</td>
        <td>61420004</td>
    </tr>
    <tr>
        <td>chrX</td>
        <td>169476592</td>
    </tr>
    <tr>
        <td>chrY</td>
        <td>91455967</td>
    </tr>
    <tr>
        <td>chrM</td>
        <td>16299</td>
    </tr>
    <tr>
        <td>GL456210.1</td>
        <td>169725</td>
    </tr>
    <tr>
        <td>GL456211.1</td>
        <td>241735</td>
    </tr>
    <tr>
        <td>GL456212.1</td>
        <td>153618</td>
    </tr>
    <tr>
        <td>GL456219.1</td>
        <td>175968</td>
    </tr>
    <tr>
        <td>GL456221.1</td>
        <td>206961</td>
    </tr>
    <tr>
        <td>GL456239.1</td>
        <td>40056</td>
    </tr>
    <tr>
        <td>GL456354.1</td>
        <td>195993</td>
    </tr>
    <tr>
        <td>GL456372.1</td>
        <td>28664</td>
    </tr>
    <tr>
        <td>GL456381.1</td>
        <td>25871</td>
    </tr>
    <tr>
        <td>GL456385.1</td>
        <td>35240</td>
    </tr>
    <tr>
        <td>JH584295.1</td>
        <td>1976</td>
    </tr>
    <tr>
        <td>JH584296.1</td>
        <td>199368</td>
    </tr>
    <tr>
        <td>JH584297.1</td>
        <td>205776</td>
    </tr>
    <tr>
        <td>JH584298.1</td>
        <td>184189</td>
    </tr>
    <tr>
        <td>JH584299.1</td>
        <td>953012</td>
    </tr>
    <tr>
        <td>JH584303.1</td>
        <td>158099</td>
    </tr>
    <tr>
        <td>JH584304.1</td>
        <td>114452</td>
    </tr>
</tbody></table>