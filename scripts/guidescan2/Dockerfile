FROM python:3.10

# Install required system packages (matching Dockerfile.processing + cellxgene-schema needs)
RUN apt-get update && \
	apt-get upgrade -y && \
	apt-get install -y \
		libvips \
		libhdf5-dev \
		python3-h5py \
		build-essential \
		pkg-config \
		coreutils \
		tabix \
		pigz \
		curl \
		wget && \
	rm -rf /var/lib/apt/lists/*

# Install guidescan-cli from GitHub releases (prebuilt binary)
# This matches the installation from Dockerfile.processing
RUN curl -fsSL https://github.com/pritykinlab/guidescan-cli/releases/download/v2.2.1/guidescan-v2.2.1-Linux.tar.gz -o /tmp/guidescan.tar.gz && \
	tar -xzf /tmp/guidescan.tar.gz -C /tmp && \
	mv /tmp/bin/guidescan /usr/local/bin/guidescan && \
	chmod +x /usr/local/bin/guidescan && \
	rm -rf /tmp/guidescan.tar.gz /tmp/bin /tmp/include

# Verify guidescan installation
RUN guidescan --version

# Install uv for fast Python package installation
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

# For lighter weight Docker images
ENV UV_NO_CACHE=1

# Copy cellxgene-schema-cli source code and requirements
COPY cellxgene_schema_cli /cellxgene_schema_cli

# Install cellxgene-schema package dependencies using uv
RUN uv pip install --system -r /cellxgene_schema_cli/requirements.txt

# Install cellxgene-schema package in editable mode
RUN uv pip install --system -e /cellxgene_schema_cli

# Verify cellxgene-schema-cli installation
RUN cellxgene-schema --help

# Set working directory for tests
WORKDIR /workspace

# Default command to show both tools are available
CMD ["sh", "-c", "echo '=== guidescan version ===' && guidescan --version && echo '\n=== cellxgene-schema-cli version ===' && cellxgene-schema --help"]

