#!/bin/bash
 
# Usage: bash download_params.sh /path/to/download/directory

# ----------------------------------------------------------

set -e # Terminate on error

# Validate parameters
if [[ $# -eq 0 ]]
  then
  echo "Error: the download directory must be provided"
  exit 1
fi
 
# Verify wget
if ! command -v wget &> /dev/null
then
    echo "Error: wget was not found; installing wget is required (sudo apt-get install wget)"
    exit 1
fi
 
# Define parameters
DOWNLOAD_DIR=$1

# Make directory
mkdir -p "${DOWNLOAD_DIR}/params"
mkdir -p "${DOWNLOAD_DIR}/params/rfdiffusion"

# Download weights
wget -P "${DOWNLOAD_DIR}/params/rfdiffusion" \
     http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt \
     http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt \
     http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt
