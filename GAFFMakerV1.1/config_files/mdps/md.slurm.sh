#!/bin/bash
#SBATCH -J AMBER_NMR_OP
#SBATCH -p multi
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --gres=gpu:1

# Decide the software version
export LD_LIBRARY_PATH=/public/software/lib/:$LD_LIBRARY_PATH
source /public/software/compiler/intel/intel-compiler-2017.5.239/bin/compilervars.sh intel64
source /public/software/profile.d/apps_gromacs_2023.2.sh

echo "Start time: $(date)"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "hostname: $(hostname)"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Job directory: $(pwd)"

# Minimization
export GMX_MAXCONSTRWARN=-1
MDRUN="gmx mdrun -nb gpu -bonded gpu -pme gpu -update gpu -gpu_id 0 "
topfile=system.top
mdpdir=../../mdps

# Energy Minimization
mkdir em 
gmx grompp -f ./mdps/em.mdp -c system.gro -r system.gro -p topol.top -o ./em/em.tpr  -n index.ndx -maxwarn 1
gmx mdrun -s ./em/em.tpr -deffnm ./em/em -ntmpi 1 -ntomp 10 -gpu_id 0 -v

# Run NVT
mkdir nvt
gmx grompp -f ./mdps/nvt.mdp -c ./em/em.gro -r ./em/em.gro -p topol.top -o ./nvt/nvt.tpr -n index.ndx -maxwarn 1
gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./nvt/nvt.tpr -deffnm ./nvt/nvt

# Run NPT
mkdir npt
gmx grompp -f ./mdps/npt.mdp -c ./nvt/nvt.gro -r ./nvt/nvt.gro -t ./nvt/nvt.cpt -p topol.top -o ./npt/npt.tpr
gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s ./npt/npt.tpr -deffnm ./npt/npt

# Run Prod MD
mkdir prod
gmx grompp -f ./mdps/md.mdp -c ./npt/npt.gro -r ./npt/npt.gro -p topol.top -o ./prod/md.tpr 
gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0 -s -deffnm ./prod/md 


echo "End time: $(date)"
