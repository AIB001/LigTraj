#!/bin/bash
# nohup bash ../../mdps/run_local.sh 2&>log & 

MDRUN="gmx mdrun -ntmpi 1 -ntomp 15 -nb gpu -bonded gpu -pme gpu -gpu_id 0"
topfile=system.top

if [ ! -f em.tpr ]; then
  gmx grompp -f ../../mdps/em.mdp -c system.gro -r system.gro -p ${topfile} -o em.tpr # -n index.ndx -maxwarn 1
fi
if [ -f em.tpr ] && [ ! -f em.gro ]; then
  gmx mdrun -s em.tpr -deffnm em -ntmpi 1 -ntomp 10 -gpu_id 0 -v
fi
if [ ! -f nvt.tpr ]; then
  gmx grompp -f ../../mdps/nvt.mdp -c em.gro -r em.gro -p ${topfile} -o nvt.tpr # -n index.ndx -maxwarn 1
fi
if [ -f nvt.tpr ] && [ ! -f nvt.gro ]; then
  $MDRUN -s nvt.tpr -deffnm nvt
fi
if [ ! -f npt.tpr ]; then
  gmx grompp -f ../../mdps/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p ${topfile} -o npt.tpr # -n index.ndx -maxwarn 1
fi
if [ -f npt.tpr ] && [ ! -f npt.gro ]; then
  $MDRUN -s npt.tpr -deffnm npt
fi
if [ ! -f md.tpr ]; then
  gmx grompp -f ../../mdps/md.mdp -c npt.gro -r npt.gro -p ${topfile} -o md.tpr # -n index.ndx -maxwarn 1
fi
if [ -f md.tpr ] && [ ! -f md.gro ]; then
  if [ -f md.cpt ]; then
    $MDRUN -deffnm md -cpi md.cpt
  else
    $MDRUN -deffnm md
  fi
fi
