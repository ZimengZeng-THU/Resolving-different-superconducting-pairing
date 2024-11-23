#!/bin/bash
#SBATCH -p paratera
#SBATCH -N 1
#SBATCH -n 64
source /public1/soft/modules/module.sh
#module load qe/6.5-intel18-zyq
module load qe/6.7
#module load qe/6.5-openmpi
#module load wannier90/intel19/3.3.0-cjj
#export PATH=/public1/home/scb0525/bin:$PATH
#srun pw.x < WTe2.relax > relax.out
#srun pw.x < scf.in > scf.out
#srun pw.x < WTe2.band > band0.out
srun ./pw.x < scf.in > scf.out
srun ./ph.x < ph.in > ph.out
#srun /public1/home/scb0525/qe/q-e-qe-7.0/bin/pw.x < scf.in > scf.out
#module load wannier90/intel19/3.3.0-cjj
#wannier90.x -pp WTe2
#srun pw2wannier90.x < WTe2.pw2wan > pw2wan.out
#wannier90.x WTe2
#yhrun pw.x < silicon.band > silicon.band.out
#srun bands.x < bands.in > bands.out
#chmod 777 plot.py
#~/.conda/envs/nlp/bin/python3.7 ./plot.py
#yhrun plotband.x < plotbands.in > plotbands.out
#yhrun wannier90.x -pp  Zn
#yhrun projwfc.x < proj.in > proj.out
#yhrun pw2wannier90.x < pw2wan.in > pw2wan.out
#yhrun wannier90.x Zn
#yhrun epw.x -npool 48  <epw.in>epw.out

#srun q2r.x <q2r.in >q2r.out
#srun matdyn.x <matdyn.in.freq >matdyn.out.freq
#yhrun matdyn.x <matdyn.in.dos>matdyn.out.dos
#yhrun lambda.x < lambda.in > lambda.out
