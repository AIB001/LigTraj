#!/bin/bash

# Usage:
# bash GAFFMaker.sh <structure_input_directory/> <simulation_file_output_director> <squence name> <expected_run_time>
# e.g. GAFFMaker.sh <structure_input_directory/> quick 12

#########################################
# Check arguments
#########################################
if [[ $# -ne 4 ]]; then
    echo "Usage: GAFFMaker.sh <structure_input_directory/> <simulation_file_output_director> <sequence name> <time>"
    exit 1
fi

structure_input_directory="${1%/}"
simulation_file_output_director="${2%/}"
seq_name="${3}"
time="${4}"

if [[ ! -d "$structure_input_directory" ]]; then
    echo "Error: structure input directory '$ligand_dir' does not exist."
    exit 1
fi

mkdir -p "${simulation_file_output_director}"

#########################################
# PDB file pre-process
#########################################

echo "Cleaning PDB file"
bash pdb_clean.sh ${structure_input_directory} "${simulation_file_output_director}"
echo "Pre-Process has done"

#########################################
# Generate ForceField File
#########################################

echo "Generating the ForceField for Simulation"
bash ligand_ffmaker.sh "${structure_input_directory}" "${simulation_file_output_director}"

#########################################
# Standardization
#########################################

echo "Standardization generated files"
bash standardization.sh "${simulation_file_output_director}"
cp -r config_files/mdps "${simulation_file_output_director}"


#########################################
# Model
#########################################

echo "GROMACS Models Building..."
bash model.sh "${simulation_file_output_director}"

#########################################
# Allocate Simulation config and submit file 
#########################################

# cp config_files/bench_submit.sh "${simulation_file_output_director}"
# cp config_files/iqb.md.slurm.sh "${simulation_file_output_director}"
# sed -i "s/SLURM_TASKNAME/${simulation_file_output_director}/g" "${simulation_file_output_director}/iqb.md.slurm.sh"
bash allocate.sh "${simulation_file_output_director}" "${seq_name}" "${time}"


#########################################
# Cleaning 
#########################################

rm sqm*
rm ANTECHAMBER*
rm PREP*
rm ATOMTYPE*
rm NEWPDB*
rm leap*
