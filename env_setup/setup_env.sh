#!/bin/bash
# take in environment name: taken from https://www.shellscript.sh/tips/getopt/

CONDA="mamba" # conda or mamba
KERNEL=false # set up jupyter kernel
INSTALL_R=true # set up R
ENV_NAME="ccc_protocols" # environment name
ENV_FILE="env_python.yml"
ENV_FILE_R="env_r.yml"
usage()
{
  echo "bash -i setup_env.sh [ -c | --conda ] [ -k | --kernel ] [ -g | --gpu ]
                             [ -p | --python-only ] [ -n | --name ENV_NAME ]"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -o ckgn: --long conda,kernel,gpu,python-only,name: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -c | --conda)  USE_CONDA=true   ; shift   ;;
    -k | --kernel) KERNEL=true      ; shift   ;;
    -g | --gpu)    GPU=true         ; shift   ;;
    -p | --python-only)   INSTALL_R=false ; shift ;;
    -n | --name)   ENV_NAME="$2" ; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done
######################################################################
# use conda or mamba? 
if [[ $USE_CONDA == true ]]
then
  CONDA="conda"
else
  MAMBA_INSTALLED="$(conda list -n base mamba | egrep "mamba")"
  if [[ ${#MAMBA_INSTALLED} == 0  ]]
  then
      echo "Mamba must be installed in your base environment, or specify option --conda"
      exit 1
  fi
fi
# GPU 
if [[ $GPU == true ]]
then
  ENV_FILE="env_python_gpu.yml"
fi

# activate environment
$CONDA env create --name "$ENV_NAME" --file $ENV_FILE
conda activate "$ENV_NAME"

# check that environment is activated
ACT_ENV="$(conda info|egrep "active environment")"
ACT_ENV=(${ACT_ENV// : / })
ACT_ENV=${ACT_ENV[2]}
if [[ "$ACT_ENV" != "$ENV_NAME" ]]; then
  echo "The environment $ENV_NAME has not been activated"
  exit 1
fi

if [[ $INSTALL_R == true ]]
then
  $CONDA env update --name "$ENV_NAME" --file $ENV_FILE_R
  conda activate "$ENV_NAME"
  Rscript r_installs.r # install the R packages not available through conda
fi

# install ipykernel and irkernel if not already installed
if [[ $KERNEL == true ]]
  then
    PKG_INSTALLED="$(conda list -n base nb_conda_kernels | egrep "nb_conda_kernels")"
    if [[ ${#PKG_INSTALLED} == 0  ]]
    then
      echo "nb_conda_kernels must be installed in your base environment, or do not specify option --kernel"
      exit 1
    fi
    PKG_INSTALLED="$(conda list -n $ENV_NAME ipykernel | egrep "ipykernel")"
    if [[ ${#PKG_INSTALLED} == 0  ]]
    then
      $CONDA install -y -c conda-forge ipykernel
    fi
    if [[ $INSTALL_R == true ]]
    then
      PKG_INSTALLED="$(conda list -n $ENV_NAME r-irkernel | egrep "r-irkernel")"
      if [[ ${#PKG_INSTALLED} == 0  ]]
      then
        $CONDA install -y -c conda-forge r-irkernel
      fi
    fi

fi
