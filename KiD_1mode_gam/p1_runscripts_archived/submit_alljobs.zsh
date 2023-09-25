#!/bin/zsh

sbatch job_condonly
sbatch job_collonly
sbatch job_sedonly
sbatch job_evaponly
sbatch job_condcoll
# sbatch job_collsed
# sbatch job_evapsed
sbatch job_condcollsed
# sbatch job_collsedevap
sbatch job_fullmic
