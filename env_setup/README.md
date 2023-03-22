# Environment Setup

This directory provides scripts to automate the creation of a conda environment for reproducible runs of the Tutorial notebooks in both python and R. 

### To create the conda environment, simply run the bash script using the following command:
```
bash -i setup_env.sh
```
*This will create the environment with both Python and R dependencies, but it does not enable a GPU.*

*The following optional flags are available as well:*
* -c | --conda: use the conda command to generate the environment, otherwise defaults to mamba
* -k | --kernel: install jupyter kernels to run notebooks in the environment
* -g | --gpu: install necessary dependencies to run a cuda-enabled GPU
* -p | --python-only: sets up only the python-associated dependencies. Use this if you do not plan using R. 
* -n | --name <ENV_NAME>: specify the environment name, otherwise defaults to "ccc_protocols"

In case of issues, .yml files are commented to include exact versions and builds used.

## Examples of other setups:
### Install Python and R dependencies enabling a GPU and kernel for jupyter notebooks (using conda):
```
bash -i setup_env.sh --conda --gpu --kernel
```

### Install only Python dependencies and enable kernel for jupyter notebooks (using conda):
```
bash -i setup_env.sh --conda --python-only --kernel
```

### Install only Python dependencies, enable a GPU  and kernel for jupyter notebooks (using conda):
```
bash -i setup_env.sh --conda --python-only --gpu --kernel
```