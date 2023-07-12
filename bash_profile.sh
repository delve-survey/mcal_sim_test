
#source your conda environment first
#Then source this file as "source bash_profile.sh"

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/oracle_instant_client/instantclient_21_1
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/despyfits/lib

export PATH=$PATH:/home/dhayaa/.local/bin/

#Hacking to get cfitsio to work :/
PATH=$PATH:/home/dhayaa/Desktop/DECADE/cfitsio-4.0.0/

#Need this to point python to desmeds master file
#export PYTHONPATH=/home/dhayaa/Desktop/DECADE/desmeds_master_erin:$PYTHONPATH
#export PYTHONPATH=/home/dhayaa/Desktop/DECADE/desmeds:$PYTHONPATH
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/mcal_sim_test:$PYTHONPATH
#export PYTHONPATH=/home/dhayaa/Desktop/DECADE/meds:$PYTHONPATH
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/WeakLensingDeblending:$PYTHONPATH
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/measure_shear/:$PYTHONPATH
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/measure_shear/metacal:$PYTHONPATH
#export PYTHONPATH=/home/dhayaa/.conda/envs/shearDM/despyfits/lib:$PYTHONPATH


export ROWE_STATS_RUN_DIR=/home/dhayaa/Desktop/DECADE/shearcat/shear_rowe_stats/
export ROWE_STATS_DIR=/scratch/midway2/dhayaa/TEMP


export PREP_DIR=/scratch/midway2/dhayaa/PREP_OUTPUTS
export TMPDIR=/scratch/midway2/dhayaa/TMP_DIR
export SHEAR_TMPDIR=/scratch/midway2/dhayaa/TMP_DIR
export MEDS_DIR=/scratch/midway2/dhayaa/MEDS_DIR
export EXP_DIR=/scratch/midway2/dhayaa/EXP_DIR
export RUN_DIR=/home/dhayaa/Desktop/DECADE/mcal_sim_test/
export MCAL_DIR=/project/chihway/dhayaa/DECADE/Tests/ 
export CATSIM_DIR=/home/dhayaa/Desktop/DECADE/WeakLensingDeblending/
export CATCOSMOS_DIR=/home/dhayaa/Desktop/DECADE/SimCatalogs/
export DESDM_CONFIG=/home/dhayaa/Desktop/DECADE/Y6DESDM
export SWARP_DIR=/home/dhayaa/Desktop/DECADE/Y6DESDM/swarp-2.40.1/
export SRCEXT_DIR=/home/dhayaa/Desktop/DECADE/Y6DESDM/sextractor-2.24.4
export LSSTSTARSIM_DIR=/scratch/midway2/dhayaa/LSSTStars/

#DESDM necessities
#export IMSUPPORT_DIR=/home/dhayaa/.conda/envs/shearDM/imsupport/
#export DESPYFITS_DIR=/home/dhayaa/.conda/envs/shearDM
#export PIXCORRECT_DIR=/home/dhayaa/Desktop/DECADE/pixcorrect

#DES-easyaccess
export DESREMOTE_RSYNC=desar2.cosmology.illinois.edu::ALLDESFiles/desarchive
export DESREMOTE_RSYNC_USER=dhayaa
export DES_RSYNC_PASSFILE=${HOME}/.desrsyncpass

export DECADEREMOTE_WGET=https://decade.ncsa.illinois.edu/deca_archive/

export EMAIL_ID=dhayaa@uchicago.edu
