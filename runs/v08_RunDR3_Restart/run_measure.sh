#!/bin/bash
#SBATCH --job-name run_measure
#SBATCH --partition=chihway
#SBATCH --account=pi-chihway
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=70:00:00
#SBATCH --output=/home/dhayaa/Desktop/DECADE/mcal_sim_test/runs/v08_RunDR3_Restart/%x.log
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --exclude=midway2-0694,midway2-0690


if [ "$USER" == "dhayaa" ]
then
    cd /home/dhayaa/Desktop/DECADE/mcal_sim_test/runs/v08_RunDR3_Restart
    module load python
    conda activate shearDM
    source /home/dhayaa/Desktop/DECADE/mcal_sim_test/bash_profile.sh
fi

python -u measure_bias.py --TileCount 1600 #--NewClassification
