This directory provides scripts to automate the creation of a conda environment for reproducible runs of the Tutorial notebooks in both python and R. 

To create the conda environment, simply run the bash script using the following command:

```
bash -i setup_env.sh
```

The following optional flags are available as well:
* -c | --conda: use the conda command to generate the environment, otherwise defaults to mamba
* -k | --kernel: install jupyter kernels to run notebooks in the environment
* -g | --gpu: install necessary dependencies to run a cuda-enabled GPU
* -n | --name <ENV_NAME>: specify the environment name, otherwise defaults to "ccc_protocols"

In case of issues, .yml files are commented to include exact versions and builds used. 