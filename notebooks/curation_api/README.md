# Discover API notebooks

These notebooks facilitate use of the [Discover API](https://api.cellxgene.cziscience.com/curation/ui/#/) which enables 
data consumers to query for public collections and datasets available in 
[CZ CELLxGENE Discover](https://cellxgene.cziscience.com/) and to download individual datasets.

## Notebook usage

Many notebooks have a few fields that require user input, such as the relevant Collection id or Dataset id.

### Python notebooks: `python/`

These notebooks are the most user-friendly; start here if you're unsure. All access token generation, url construction, and response 
handling is abstracted away from the user. Unitary interactions such as getting a Dataset, getting a Collection, etc., 
are each handled by single function calls.

### "Raw" Python notebooks: `python_raw/`

These notebooks show the construction of http requests at a more granular level and may be useful for technical
users.

### R notebooks: `R/`

These notebooks are the same as the "raw" Python notebooks mentioned immediately previous, but in R.
