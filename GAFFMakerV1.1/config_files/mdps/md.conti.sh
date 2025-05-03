#!/bin/bash
#SBATCH -J run4-conti
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

mdpdir=../../mdps
# Minimization
export GMX_MAXCONSTRWARN=-1
MDRUN="gmx mdrun -nb gpu -bonded gpu -pme gpu -update gpu -gpu_id 0 "

if [ ! -f md2.tpr ]; then
  gmx grompp -f ${mdpdir}/md.mdp -c md.gro -r md.gro -p topol.top -o md2.tpr -n index.ndx -maxwarn 1
fi
if [ -f md2.tpr ] && [ ! -f md2.gro ]; then
  if [ -f md2.cpt ]; then
    $MDRUN -deffnm md2 -cpi md2.cpt
  else
    $MDRUN -deffnm md2
  fi
fi

echo "End time: $(date)"
