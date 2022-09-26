# STATE-BATCH
Python wrapper for high-throughput calculations using STATE code.

**Dependencies:**\
Below are dependencies that should be installed before running STATE-BATCH.
- ASE-STATE interface 
  - `git clone -b fix-bug-output-parser git@github.com:uedar/ase-state-interface.git`
  - `cd ase-state-interface`
  - `pip install .`
- pandas 
  - `pip install pandas`
- yaml
  - `pip install PyYaml`

**How to run:**\
Prepare three files, namely: input.yaml, a csv file, and run.py in one directory.\
Run the python code: `python run.py`

**Details:**\
Those files are explained below:

- **input.yaml**\
This yaml file contains all the settings necessary for the high throughput calculations.\
There are three important components need to be specified: `comp_spec`, `system_spec`, and `dft_spec`.\
Please see examples and the comments within for more details.

- **csv file**\
The csv file contains the list of the system to be run in high-throughput fashion.

- **run.py**\
This python code is used to run STATE-BATCH.
