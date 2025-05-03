#!/bin/bash
#SBATCH -J SIRPa-pep-oplsaam-restr
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

if [ ! -f em.tpr ]; then
  gmx grompp -f ${mdpdir}/em.mdp -c system.gro -r system.gro -p ${topfile} -o em.tpr # -n index.ndx -maxwarn 1
fi
if [ -f em.tpr ] && [ ! -f em.gro ]; then
  gmx mdrun -s em.tpr -deffnm em -ntmpi 1 -ntomp 10 -gpu_id 0 -v
fi
if [ ! -f nvt.tpr ]; then
  gmx grompp -f ${mdpdir}/nvt.mdp -c em.gro -r em.gro -p ${topfile} -o nvt.tpr # -n index.ndx -maxwarn 1
fi
if [ -f nvt.tpr ] && [ ! -f nvt.gro ]; then
  $MDRUN -s nvt.tpr -deffnm nvt
fi
if [ ! -f npt.tpr ]; then
  gmx grompp -f ${mdpdir}/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ${topfile} -o npt.tpr # -n index.ndx -maxwarn 1
fi
if [ -f npt.tpr ] && [ ! -f npt.gro ]; then
  $MDRUN -s npt.tpr -deffnm npt
fi
if [ ! -f md_restr.tpr ]; then
  gmx grompp -f ${mdpdir}/md_restr.mdp -c npt.gro -r npt.gro -p ${topfile} -o md_restr.tpr # -n index.ndx -maxwarn 1
fi
if [ -f md_restr.tpr ] && [ ! -f md_restr.gro ]; then
  if [ -f md_restr.cpt ]; then
    $MDRUN -deffnm md_restr -cpi md_restr.cpt
  else
    $MDRUN -deffnm md_restr
  fi
fi
if [ ! -f md.tpr ]; then
  gmx grompp -f ${mdpdir}/md.mdp -c md_restr.gro -r md_restr.gro -p ${topfile} -o md.tpr # -n index.ndx -maxwarn 1
fi
if [ -f md.tpr ] && [ ! -f md.gro ]; then
  if [ -f md.cpt ]; then
    $MDRUN -deffnm md -cpi md.cpt
  else
    $MDRUN -deffnm md
  fi
fi

echo "End time: $(date)"
