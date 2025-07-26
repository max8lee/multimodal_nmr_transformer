#!/bin/bash

# Create a data directory if it doesn't exist
mkdir -p data

# Download the data from Zenodo
echo "Downloading spectra data from Zenodo..."
curl -L "https://zenodo.org/records/14770232/files/multimodal_spectroscopic_dataset.zip" -o data/multimodal_spectroscopic_dataset.zip

# Unzip the file
echo "Extracting data..."
UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE unzip -q data/multimodal_spectroscopic_dataset.zip -d data/

# Remove the zip file to save space (optional)
# echo "Cleaning up..."
# rm data/multimodal_spectroscopic_dataset.zip

echo "Download complete! Data is available in the data directory."