---
name: Add species
about: Editor's template for adding new species
title: Draft <species>
labels: drafting, multispecies discovery, schema
assignees: brianraymor

---

## Pending Issues

1. Waiting on sscrdv to be submitted to OLS for use in references
1. [FAANG](http://www.faang.org/) is the Functional Annotation of ANimal Genomes project. _We are working to understand the genotype to phenotype link in domesticated animals._ Per their [Ontology Improver](https://data.faang.org/ontology?sortTerm=key&sortDirection=asc), *Dv terms are not referenced. Both UBERON and CL are in use. Their [schema](https://github.com/FAANG/dcc-metadata/blob/9e7c1b5304fc57a724d197384e83243562bebbf4/json_schema/type/samples/faang_samples_specimen.metadata_rules.json#L154):

```
"name": "developmental stage",
"description": "Ontology for Developmental stage, UBERON is preferred to EFO.",
```


## Design

This draft design reflects additions to corresponding sections in [schema 5.2.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.2.0/schema.md). Reviewers are expected to be familiar with the CELLxGENE schema.

**Editorial Notes** that are inlined in the design below will not be surfaced in the schema. 

---

### Required Ontologies


| Ontology | OBO Prefix | Release | Download |
|:--|:--|:--|:--|
| [Unavailable](https://github.com/OBOFoundry/OBOFoundry.github.io/tree/master/ontology) | SscrDv | [Releases](https://github.com/obophenotype/developmental-stage-ontologies/releases) | TBD |
|||||


#### Editorial Notes

This ontology is under active development. CELLxGENE pins ontology releases in each version of the schema. A specific release of the ontology above must be selected in the future.


---

### Required Gene Annotations

| Organism | Source | Required version | Download |
|:--|:--|:--|:--|
| <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9823"><code>"NCBITaxon:9823"</code></a><br>for <i>Sus scrofa domesticus</i>  | [ENSEMBL (Sus scrofa domesticus)] | Sscrofa11.1 (GCA_000003025.6) | [Sus_scrofa.Sscrofa11.1.113.gtf] |


[ENSEMBL (Sus scrofa domesticus)]: https://useast.ensembl.org/Sus_scrofa/Info/Index
[Sus_scrofa.Sscrofa11.1.113.gtf]: https://ftp.ensembl.org/pub/release-113/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.113.gtf.gz

#### Editorial Notes


---

## `obs` (Cell Metadata)

### cell_type_ontology_term_id

No schema changes are required. 

#### Editorial Notes

---

### development_stage_ontology_term_id

<table><tbody>
  <tr>
    <th>Key</th>
    <td>development_stage_ontology_term_id</td>
  </tr>
  <tr>
    <th>Annotator</th>
    <td>Curator MUST annotate.</td>
  </tr>
  <tr>
    <th>Value</th>
    <td>
      categorical with <code>str</code> categories. If unavailable, this MUST be <code>"unknown"</code>.<br><br>
      If <code>organism_ontolology_term_id</code> is <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9823"><code>"NCBITaxon:9823"</code></a> for <i>Sus scrofa domesticus</i>, this MUST be the most accurate descendant of <code>SscrDv:0000000</code> for <i>life cycle stage<i>.
    </td>
  </tr>
</tbody></table>
<br>

#### Editorial Notes

This may be outdated, but [potential recommendations](https://github.com/obophenotype/developmental-stage-ontologies/blob/master/external/bgee/report.md#sus-scrofa):

```
UBERON:0000104 life cycle
    UBERON:0000068 embryo stage
        UBERON:0000106 zygote stage
        UBERON:0000107 cleavage stage
            UBERON:0007232 2 cell stage
            UBERON:0007233 4 cell stage
            UBERON:0007236 8 cell stage
        UBERON:0000108 blastula stage
        UBERON:0000109 gastrula stage
        UBERON:0000110 neurula stage
        UBERON:0000111 organogenesis stage
            SscrDv:0000081 ridge limb stage (pig)
            SscrDv:0000082 bud limb stage (pig)
            SscrDv:0000083 paddle limb stage (pig)
        UBERON:0007220 late embryonic stage
    UBERON:0000092 post-embryonic stage
        UBERON:0000066 fully formed stage
            UBERON:0000112 sexually immature stage
                UBERON:0018685 nursing stage
                    UBERON:0007221 neonate stage
                        SscrDv:0000072 0-day-old stage (pig)
                        SscrDv:0000073 1-day-old stage (pig)
                        SscrDv:0000074 2-day-old stage (pig)
                        SscrDv:0000075 3-day-old stage (pig)
                        SscrDv:0000076 4-day-old stage (pig)
                        SscrDv:0000077 5-day-old stage (pig)
                        SscrDv:0000078 6-day-old stage (pig)
                    UBERON:0034920 infant stage
                        SscrDv:0000010 1-week-old stage (pig)
                        SscrDv:0000011 2-week-old stage (pig)
                        SscrDv:0000012 3-week-old stage (pig)
                            SscrDv:0000018 21-day-old stage (pig)
                            SscrDv:0000019 22-day-old stage (pig)
                            SscrDv:0000020 23-day-old stage (pig)
                            SscrDv:0000021 24-day-old stage (pig)
                            SscrDv:0000022 25-day-old stage (pig)
                            SscrDv:0000023 26-day-old stage (pig)
                            SscrDv:0000024 27-day-old stage (pig)
                        SscrDv:0000013 4-week-old stage (pig)
                            SscrDv:0000025 28-day-old stage (pig)
                            SscrDv:0000026 29-day-old stage (pig)
                            SscrDv:0000027 30-day-old stage (pig)
                            SscrDv:0000028 31-day-old stage (pig)
                            SscrDv:0000029 32-day-old stage (pig)
                            SscrDv:0000030 33-day-old stage (pig)
                            SscrDv:0000031 34-day-old stage (pig)
                        SscrDv:0000014 5-week-old stage (pig)
                            SscrDv:0000032 35-day-old stage (pig)
                            SscrDv:0000033 36-day-old stage (pig)
                            SscrDv:0000034 37-day-old stage (pig)
                            SscrDv:0000035 38-day-old stage (pig)
                            SscrDv:0000036 39-day-old stage (pig)
                            SscrDv:0000037 40-day-old stage (pig)
                            SscrDv:0000038 41-day-old stage (pig)
                        SscrDv:0000015 6-week-old stage (pig)
                        SscrDv:0000016 7-week-old stage (pig)
                UBERON:0034919 juvenile stage
                    SscrDv:0000039 2-month-old stage (pig)
                        SscrDv:0000017 8-week-old stage (pig)
                        SscrDv:0000040 9-week-old stage (pig)
                        SscrDv:0000041 10-week-old stage (pig)
                        SscrDv:0000042 11-week-old stage (pig)
                    SscrDv:0000043 3-month-old stage (pig)
                        SscrDv:0000044 12-week-old stage (pig)
                        SscrDv:0000045 13-week-old stage (pig)
                        SscrDv:0000046 14-week-old stage (pig)
                        SscrDv:0000047 15-week-old stage (pig)
                    SscrDv:0000048 4-month-old stage (pig)
                        SscrDv:0000049 16-week-old stage (pig)
                        SscrDv:0000050 17-week-old stage (pig)
                        SscrDv:0000051 18-week-old stage (pig)
                        SscrDv:0000052 19-week-old stage (pig)
                        SscrDv:0000053 20-week-old stage (pig)
                    SscrDv:0000054 5-month-old stage (pig)
                        SscrDv:0000055 21-week-old stage (pig)
                        SscrDv:0000056 22-week-old stage (pig)
                        SscrDv:0000057 23-week-old stage (pig)
                        SscrDv:0000058 24-week-old stage (pig)
                    SscrDv:0000059 6-month-old stage (pig)
                    SscrDv:0000060 7-month-old stage (pig)
                    SscrDv:0000061 8-month-old stage (pig)
                    SscrDv:0000062 9-month-old stage (pig)
                    SscrDv:0000063 10-month-old stage (pig)
            UBERON:0000113 post-juvenile
                UBERON:0018241 prime adult stage
                    SscrDv:0000064 11-month-old stage (pig)
                    SscrDv:0000065 1-year-old stage (pig)
                    SscrDv:0000066 2-year-old stage (pig)
                    SscrDv:0000067 3-year-old stage (pig)
                    SscrDv:0000068 4-year-old stage (pig)
                    SscrDv:0000069 5-year-old stage (pig)
                    SscrDv:0000070 6-year-old stage (pig)
                    SscrDv:0000071 7-year-old stage (pig)
                UBERON:0007222 late adult stage
```

---

### disease_ontology_term_id

No schema changes are required.

#### Editorial Notes

---

### organism_ontolology_term_id

<code>organism_ontolology_term_id</code> is <a href="https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon/classes?obo_id=NCBITaxon%3A9823"><code>"NCBITaxon:9823"</code></a> for <i>Sus scrofa domesticus</i> 

---

### sex_ontology_term_id

No schema changes are required.

#### Editorial Notes

---

### tissue_ontology_term_id

No schema changes are required.


#### Editorial Notes

---

## Reference


[BGEE](https://www.bgee.org/species/9823)
