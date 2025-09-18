#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# ===================================================================
# Configuration Variables
# ===================================================================
# Define the base directory for all project-related files.
# This makes the script portable. All other paths are relative to this one.
BASE_DIR=".."

# Define sub-directories based on the base path.
DATA_DIR="${BASE_DIR}/Data"
RESULT_DIR="${BASE_DIR}/result"
CODE_DIR="${BASE_DIR}/code"

# Define variables for specific datasets and results to keep paths clean.
SLIDE1_RAW_DATA_PATH="${DATA_DIR}/zfnkx25-14ss-slide1-020924_VMSC02901/region_"
SLIDE2_RAW_DATA_PATH="${DATA_DIR}/zfnkx25-14ss-slide2-021224_VMSC02901/region_"

# Define result paths for the requested segmentation method.
SLIDE1_CELLPOSE_RESULT_PATH="${RESULT_DIR}/cell_seg_DAPI_mouse_rabbit/slide1region"
SLIDE2_CELLPOSE_RESULT_PATH="${RESULT_DIR}/cell_seg_DAPI_mouse_rabbit/slide2region"

# Define path to the segmentation algorithm.
CELLPOSE_ALGORITHM_PATH="${CODE_DIR}/example_analysis_algorithm/cellpose_default_1_ZLevel_custom01.json"

# List the regions to be processed.
SLIDE1_REGIONS=(0 2 3)
SLIDE2_REGIONS=(0 2 4 5)

# ===================================================================
# Conda Environment Setup
# ===================================================================
# Activate the correct conda environment for the analysis.
source ~/.bashrc
conda activate vpt_cellpose2

date

# ===================================================================
# Main Analysis Pipeline (Slide 1 - Cellpose)
# ===================================================================
echo "Starting cell segmentation pipeline for Slide 1..."

# Loop through each region for Slide 1.
for i in "${SLIDE1_REGIONS[@]}"; do
    echo "Processing Slide 1, Region: ${i}"

    # Create the output directory for the current region if it doesn't exist.
    mkdir -p "${SLIDE1_CELLPOSE_RESULT_PATH}${i}"

    # Step 1: Identify Cell Boundaries
    vpt --verbose --processes 4 run-segmentation \
        --segmentation-algorithm "${CELLPOSE_ALGORITHM_PATH}" \
        --input-images="${SLIDE1_RAW_DATA_PATH}${i}/images/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" \
        --input-micron-to-mosaic "${SLIDE1_RAW_DATA_PATH}${i}/images/micron_to_mosaic_pixel_transform.csv" \
        --output-path "${SLIDE1_CELLPOSE_RESULT_PATH}${i}" \
        --tile-size 2400 \
        --tile-overlap 200

    # Step 2: Partition Transcripts into Cells
    vpt --verbose partition-transcripts \
        --input-boundaries "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cellpose_micron_space.parquet" \
        --input-transcripts "${SLIDE1_RAW_DATA_PATH}${i}/detected_transcripts.csv" \
        --output-entity-by-gene "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cell_by_gene.csv" \
        --output-transcripts "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/detected_transcripts.csv"

    # Step 3: Calculate Cell Metadata and Sum Signals
    vpt --verbose derive-entity-metadata \
        --input-boundaries "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cellpose_micron_space.parquet" \
        --input-entity-by-gene "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cell_by_gene.csv" \
        --output-metadata "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cell_metadata.csv"

    vpt --verbose sum-signals \
        --input-images="${SLIDE1_RAW_DATA_PATH}${i}/images/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" \
        --input-boundaries "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cellpose_micron_space.parquet" \
        --input-micron-to-mosaic "${SLIDE1_RAW_DATA_PATH}${i}/images/micron_to_mosaic_pixel_transform.csv" \
        --output-csv "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/sum_signals.csv"

    # Step 4: Update the .vzg File
    vpt --verbose --processes 2 update-vzg \
        --input-vzg "${SLIDE1_RAW_DATA_PATH}${i}/202402091240_hakan-zfnkx25-14ss-slide1020924_VMSC02901_region_${i}.vzg" \
        --input-boundaries "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cellpose_micron_space.parquet" \
        --input-entity-by-gene "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cell_by_gene.csv" \
        --input-metadata "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/cell_metadata.csv" \
        --output-vzg "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/slide1_region_${i}.vzg"

    # Copy the micron_to_mosaic file to the result directory.
    mkdir -p "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/images"
    cp "${SLIDE1_RAW_DATA_PATH}${i}/images/micron_to_mosaic_pixel_transform.csv" "${SLIDE1_CELLPOSE_RESULT_PATH}${i}/images/"

    echo "Finished processing Slide 1, Region: ${i}"
    echo "---------------------------------------------------"
done

# ===================================================================
# Main Analysis Pipeline (Slide 2 - Cellpose)
# ===================================================================
echo "Starting cell segmentation pipeline for Slide 2..."

# Loop through each region for Slide 2.
for i in "${SLIDE2_REGIONS[@]}"; do
    echo "Processing Slide 2, Region: ${i}"

    # Create the output directory for the current region if it doesn't exist.
    mkdir -p "${SLIDE2_CELLPOSE_RESULT_PATH}${i}"

    # Step 1: Identify Cell Boundaries
    vpt --verbose --processes 4 run-segmentation \
        --segmentation-algorithm "${CELLPOSE_ALGORITHM_PATH}" \
        --input-images="${SLIDE2_RAW_DATA_PATH}${i}/images/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" \
        --input-micron-to-mosaic "${SLIDE2_RAW_DATA_PATH}${i}/images/micron_to_mosaic_pixel_transform.csv" \
        --output-path "${SLIDE2_CELLPOSE_RESULT_PATH}${i}" \
        --tile-size 2400 \
        --tile-overlap 200

    # Step 2: Partition Transcripts into Cells
    vpt --verbose partition-transcripts \
        --input-boundaries "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cellpose_micron_space.parquet" \
        --input-transcripts "${SLIDE2_RAW_DATA_PATH}${i}/detected_transcripts.csv" \
        --output-entity-by-gene "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cell_by_gene.csv" \
        --output-transcripts "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/detected_transcripts.csv"

    # Step 3: Calculate Cell Metadata and Sum Signals
    vpt --verbose derive-entity-metadata \
        --input-boundaries "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cellpose_micron_space.parquet" \
        --input-entity-by-gene "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cell_by_gene.csv" \
        --output-metadata "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cell_metadata.csv"

    vpt --verbose sum-signals \
        --input-images="${SLIDE2_RAW_DATA_PATH}${i}/images/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" \
        --input-boundaries "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cellpose_micron_space.parquet" \
        --input-micron-to-mosaic "${SLIDE2_RAW_DATA_PATH}${i}/images/micron_to_mosaic_pixel_transform.csv" \
        --output-csv "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/sum_signals.csv"

    # Step 4: Update the .vzg File
    vpt --verbose --processes 2 update-vzg \
        --input-vzg "${SLIDE2_RAW_DATA_PATH}${i}/202402121215_hakan-zfnkx25-slide2-021224_VMSC02901_region_${i}.vzg" \
        --input-boundaries "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cellpose_micron_space.parquet" \
        --input-entity-by-gene "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cell_by_gene.csv" \
        --input-metadata "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/cell_metadata.csv" \
        --output-vzg "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/slide2_region_${i}.vzg"

    # Copy the micron_to_mosaic file to the result directory.
    mkdir -p "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/images"
    cp "${SLIDE2_RAW_DATA_PATH}${i}/images/micron_to_mosaic_pixel_transform.csv" "${SLIDE2_CELLPOSE_RESULT_PATH}${i}/images/"

    echo "Finished processing Slide 2, Region: ${i}"
    echo "---------------------------------------------------"
done

echo "Pipeline complete."
date
