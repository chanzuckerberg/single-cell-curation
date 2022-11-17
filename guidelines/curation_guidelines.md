# Guidelines for CELLxGENE Curators

Identifying Data to Curate
----------------
Once a curator identifies a study for CELLxGENE submission, they should also identify any secondary DOIs, the publication DOIs corresponding to any data that the primary study includes. The curator should then check the [Dataset sheet](https://docs.google.com/spreadsheets/d/1ax9b5sxmxSJgrjncXG5WGilgIGKm2EEWmILpoN6pLzY/edit?usp=sharing) for the primary study DOI to see if it has already been submitted or is currently being curated. The curator should also check the [reused data tab](https://docs.google.com/spreadsheets/d/1ax9b5sxmxSJgrjncXG5WGilgIGKm2EEWmILpoN6pLzY/edit?usp=sharing) for the primary study DOI and both tabs for each secondary DOI to see if any data from the study have been included in other CELLxGENE submissions. Data that are already included can likely still be submitted, but the data redundancy should be discussed with the Lattice team as it will have implications on how to curate [`is_primary_data`](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md#is_primary_data). The curator should email the [Lattice team](mailto:lattice-info@lists.stanford.edu) with the DOI, as well as the secondary DOIs, so that the sheets can be updated and any data redundancies can be addressed.

Publishing Collections
----------------
Once the Collection has been created, all Datasets have been uploaded, and both the curator and contributor consider the submission finalized, the curator should email the [Lattice team](mailto:lattice-info@lists.stanford.edu) with a link to the private Collection. We ask curators to send one email per Collection as a way to easily track discussions and requested updates that may come up.

The Lattice team will review the data and may ask clarifying questions or suggest revisions before approving the Collection. Once the Lattice team has given approval and the contributor is ready to Publish, the curator may Publish the Collection. If any updates are made after the Lattice approval, then the curator should contact Lattice again for review and provide details on the scope of the updates.

Upon Publishing, the CELLxGENE communication team asks that the following message be sent to the contributor (after filling in the `Collection URL`):

>The CELLxGENE team would like to maximize the impact of your research by spreading the word that it’s now available to explore online using CZ CELLxGENE Discover. Help us by sharing on your social channels or with your institutions. A few tips:
>- Share the direct link to your collection: `Collection URL`?utm_campaign=partner&utm_source=publisher
>- Tag your post with the hashtag #CZCELLxGENE. Our social team monitors this hashtag and we may retweet your post on CZI’s channels so it reaches an even bigger audience.

Publishing Revisions
----------------
When making changes to a Dataset, the curator should Update the relevant Dataset, rather than Adding a new one. Update will retain the Dataset id and thus, the visualization URL, when Published.

If there are any questions as to whether the changes constitute a Dataset Revision or a new Dataset, please contact the [Lattice team](mailto:lattice-info@lists.stanford.edu) or raise this on the biweekly curator call for discussion.

As with Publishing new Collections, the curator should email the [Lattice team](mailto:lattice-info@lists.stanford.edu) with a link to the Revision for review, and may Publish once approved.

Deleting Published Datasets
----------------
Removing Data that have been publicly available should be avoided, but it is possible.
Please contact the [Lattice team](mailto:lattice-info@lists.stanford.edu) with any Deletion requests along with a thorough justification for the removal.
