FROM ubuntu:latest
LABEL authors="trentsmith"

# For lighter weight Docker images
ENV PIP_NO_CACHE_DIR=1


RUN apt-get update \
    && apt-get upgrade -y \
    && apt install software-properties-common -y

RUN add-apt-repository ppa:deadsnakes/ppa && \
	apt-get update && \
    apt-get install -y python3.10 python3.10-dev python3.10-venv build-essential libvips libhdf5-dev python3-h5py wget pkg-config && \
    wget https://bootstrap.pypa.io/get-pip.py && \
    python3.10 get-pip.py \
	&& apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Activate virtual environment for subsequent commands
RUN python3.10 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install test dependencies
RUN pip install memray memory-profiler

ADD ../../cellxgene_schema_cli /cellxgene_schema_cli
RUN	pip install -e /cellxgene_schema_cli

ADD ./scripts/memtest /memtest
