# HCA Curation Fork

Clever Canary's fork of [cellxgene-schema](https://github.com/chanzuckerberg/single-cell-curation) with Human Cell Atlas (HCA) validation extensions.

## HCA Schema Validator Package

**📦 `hca-schema-validator`** - Python package extending cellxgene-schema with HCA-specific validation rules.

**📖 Complete documentation:** [hca_schema_validator/README.md](hca_schema_validator/README.md)

For installation, usage, development, and testing instructions, see the package README.

## About This Fork

This fork contains:

- **`hca_schema_validator/`** - HCA-specific validation package (published to PyPI)
- **`cellxgene_schema_cli/`** - Upstream cellxgene-schema (synced from CZI)

## Syncing from Upstream

We maintain two branches:

- `main` - Our HCA customizations
- `upstream-main` - Tracks CZI's main branch

To sync updates from upstream:

```bash
git checkout upstream-main
git pull upstream main
git checkout main
git merge upstream-main
```

## Contributing

See individual package READMEs:

- [HCA Schema Validator Contributing](hca_schema_validator/README.md)
- [CellxGene Schema Contributing](README.md)

## License

MIT (same as upstream)
