#!/bin/bash

# Parse command-line arguments
while getopts ":c:" opt; do
  case $opt in
    c) config_file="$OPTARG"
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
    ;;
  esac
done

dataset=$(jq -r '.dataset' "$config_file")

cd /cbica/projects/luo_wm_dev/input/${dataset}/datalad_qsiprep
missing=$(git annex find --not --in here)
datalad get $missing 
 