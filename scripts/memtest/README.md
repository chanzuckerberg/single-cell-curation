# Memory Test
The purpose of this script to profile the memory usage of the cellxgene-schema-cli validation and label writing process.

There are two library that can be used to profile. Using both is recommended to get a more complete picture of the memory usage. The tests are run in a docker container to constrain the memory availible to the process. This is done to simulate the possible memory constraints of the ingest pipeline.

# How to run the memory test
To use either profiler you must first follow these steps:
1. Get an unlabled dataset
2. Copy the dataset to the ./data folder
3. Update H5AD_FILE in run.py to point to the dataset.
4. Install memory profiling tools in your local environment 
```bash
    pip install -r ./requirements.txt
```
5. build the docker container
```bash
docker compose --project-directory ../../ build memtest
```

## [memray](https://github.com/bloomberg/memray) 
To run the test with memray follow these steps:
1. Run the docker container.
```bash
docker compose --project-directory ../../ run --rm memtest memray run /memtest/run.py
```
2. Check the results
``` bash
memray flamegraph ${results_file}
```

## [memory-profiler](https://pypi.org/project/memory-profiler/)
To run the test with memory-profiler follow these steps:

1. run the docker container.
```bash
docker compose --project-directory ../../ run --rm memtest mprof run /memtest/run.py
```
2. Check the results
``` bash
mprof plot
```
