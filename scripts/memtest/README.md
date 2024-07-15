# How to run the memory test
This tests uses [memray](https://github.com/bloomberg/memray) to profile the memory usage of a python program. To run the test follow these steps:
1. Get an unlabled dataset
2. install the requirements
3. Update H5AD_FILE in run.py to point to the dataset.
4. run the test with the following command
```bash
python3 -m memray run run.py
```
5. Check the results
``` bash
memray flamegraph ${results_file}
```
