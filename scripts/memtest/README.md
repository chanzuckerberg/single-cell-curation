# How to run the memory test
This tests uses [memray](https://github.com/bloomberg/memray) to profile the memory usage of a python program. To run the test follow these steps:
1. Get an unlabled dataset
2. copy the dataset to the ./data folder
3. Update H5AD_FILE in run.py to point to the dataset.
4. run the docker container.
```bash
docker-compose -f ./single-cell-curation/scripts/memtest/docker-compose.yml -p memtest up -d processing
```
A docker container is used to constrain the memory usage of the program.
5. Check the results
``` bash
memray flamegraph ${results_file}
```
