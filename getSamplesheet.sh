#!/bin/bash
# Download the CSV file from a URL
wget "https://docs.google.com/spreadsheets/d/1jDX2pudze8Bk_DxtwtOrKYF-oBXHCzjPtykY5JY-6bI/export?format=csv" -O "samplesheet.csv"
current_dir=$(pwd)

# Convert the CSV file to JSON using a Python script
python $current_dir/project/xsvato01/mareckova_CXCR_k8s/src/CsvToJson.py $current_dir/samplesheet.csv $current_dir/samplesheet.json
rm samplesheet.csv