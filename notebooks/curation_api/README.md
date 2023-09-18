# Discover API notebooks

The Discover API notebooks are pre-baked scripts to facilitate use of the 
[Discover API](https://api.cellxgene.cziscience.com/curation/ui/#/) (formerly known as the Curation API), which itself supports 
curation flow interactions with the [CZ CELLxGENE](https://cellxgene.cziscience.com/) tool.

## Notebook usage

Many notebooks have a few fields that require user input, such as the relevant Collection id or Dataset id. For 
privileged user interactions (mutations), be sure to provide the path to your API key where specified in the notebook.

### Python notebooks: `python/`

These notebooks are the most user-friendly; start here if you're unsure. All access token generation, url construction, and response 
handling is abstracted away from the user. Unitary interactions such as getting a Dataset, creating a Collection, etc., 
are each handled by single function calls.

### "Raw" Python notebooks: `python_raw/`

These notebooks show the construction of http requests at a more granular level, and as such may be useful for technical
users.

### R notebooks: `R/`

These notebooks are the same as the "raw" Python notebooks mentioned immediately previous, but in R.
