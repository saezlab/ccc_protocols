#!/bin/bash
# take in environment name: taken from https://www.shellscript.sh/tips/getopt/

NAME="ccc_protocols"
usage()
{
  echo "Usage: bash -i setup_env.sh [ -n | --name ENV_NAME ]"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -o n: --long name: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -n | --name) NAME="$2" ; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

conda env create --name "$NAME" --file env.yml
conda activate "$NAME"

# check that environment is activated
ACT_ENV="$(conda info|egrep "active environment")"
ACT_ENV=(${ACT_ENV// : / })
ACT_ENV=${ACT_ENV[2]}
if [[ "$ACT_ENV" != "$NAME" ]]; then
  echo "The environment $NAME has not been activated"
fi

Rscript setup_env.r # install the R packages not available through conda