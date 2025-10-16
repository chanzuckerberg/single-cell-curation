# HCA Schema Validator

HCA-specific extensions for cellxgene schema validation.

## Installation

### From PyPI (Recommended)

```bash
pip install hca-schema-validator
```

### From Source (Development)

### Prerequisites

Install [Poetry](https://python-poetry.org/docs/#installation) if you haven't already:
```bash
curl -sSL https://install.python-poetry.org | python3 -
```

### Install

```bash
# Install cellxgene-schema dependency (uses pip)
cd ../cellxgene_schema_cli
pip install -e .

# Install hca-schema-validator (uses Poetry)
cd ../hca_schema_validator
poetry install
```

## Usage

```python
from hca_schema_validator import HCAValidator

# Create validator instance
validator = HCAValidator()

# Validate an h5ad file
is_valid = validator.validate_adata("path/to/file.h5ad")

# Check results
if is_valid:
    print("✅ Validation passed!")
else:
    print("❌ Validation failed:")
    for error in validator.errors:
        print(f"  - {error}")
```

## Development Status

**Current Version: 0.1.0** - Minimal passthrough implementation

Currently a passthrough wrapper around cellxgene-schema Validator.
HCA-specific validation rules will be added incrementally.

## Testing

```bash
cd hca_schema_validator
poetry run pytest tests/
```

## Project Structure

```
hca_schema_validator/
├── src/
│   └── hca_schema_validator/
│       ├── __init__.py       # Package exports
│       └── validator.py      # HCAValidator class
├── tests/
│   └── test_validator.py # Unit tests
├── pyproject.toml        # Poetry configuration & dependencies
└── README.md            # This file
```

## License

MIT
