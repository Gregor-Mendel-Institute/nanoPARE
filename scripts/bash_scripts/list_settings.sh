#!/bin/bash

echo "
_____________________________
SETTINGS
-----------------------------
Basic configuration:
  root_dir="$root_dir"
  bash_dir="$bash_dir"
  python_dir="$python_dir"
  resource_dir="$resource_dir"
  temp_dir="$temp_dir"
  log_dir="$log_dir"
  results_dir="$results_dir"

General settings:
  GENOME_FASTA="$GENOME_FASTA"
  ANNOTATION_GFF="$ANNOTATION_GFF"
  LMOD="$LMOD"
  RAM="$RAM"
  CPUS="$CPUS"
  
EndMap settings:
  LINE_NUMBER="$JOB_NUMBER"
  ICOMP="$ICOMP"

EndGraph settings:
  SAMPLE_NAME="$SAMPLE_NAME"
  KERNEL="$KERNEL"

EndClass settings:
  SAMPLE_TYPE="$SAMPLE_TYPE"
  UUG="$UUG"

EndMask settings:
  SAMPLE_TYPE="$SAMPLE_TYPE"
  MASK_SOURCE="$MASK_SOURCE"
_____________________________

"
