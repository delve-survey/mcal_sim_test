
#source your conda environment first
#Then source this file as "source bash_profile.sh"

export PATH=$PATH:/home/dhayaa/.local/bin/

#Hacking to get cfitsio to work :/
PATH=$PATH:/home/dhayaa/Desktop/DECADE/cfitsio-4.0.0/

#Need this to point python to desmeds master file
#export PYTHONPATH=/home/dhayaa/Desktop/DECADE/desmeds:$PYTHONPATH
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/mcal_sim_test:$PYTHONPATH
#export PYTHONPATH=/home/dhayaa/Desktop/DECADE/meds:$PYTHONPATH
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/WeakLensingDeblending:$PYTHONPATH
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/measure_shear/:$PYTHONPATH
export PYTHONPATH=/home/dhayaa/Desktop/DECADE/measure_shear/metacal:$PYTHONPATH


export ROWE_STATS_RUN_DIR=/home/dhayaa/Desktop/DECADE/shearcat/shear_rowe_stats/
export ROWE_STATS_DIR=/scratch/midway2/dhayaa/TEMP

#PATH for files
export EBV_PATH=/project/chihway/dhayaa/DECADE/Imsim_Inputs/ebv_sfd98_fullres_nside_4096_ring_equatorial.fits
export CATSIM_PATH=/project/chihway/dhayaa/DECADE/Imsim_Inputs/OneDegSq.fits
export CATCOSMOS_PATH=/project/chihway/dhayaa/DECADE/Imsim_Inputs/input_cosmos_v4.fits
export CATDESDF_PATH=/project/chihway/dhayaa/DECADE/Imsim_Inputs/DESY3_Deepfields_V1BalrogCuts.fits

export PREP_DIR=/scratch/midway2/dhayaa/PREP_OUTPUTS
export TMPDIR=/scratch/midway2/dhayaa/TMP_DIR
export SHEAR_TMPDIR=/scratch/midway2/dhayaa/TMP_DIR
export MEDS_DIR=/scratch/midway2/dhayaa/MEDS_DIR
export EXP_DIR=/scratch/midway2/dhayaa/EXP_DIR
export RUN_DIR=/home/dhayaa/Desktop/DECADE/mcal_sim_test/
export MCAL_DIR=/project/chihway/dhayaa/DECADE/Tests/ 
export DESDM_CONFIG=/home/dhayaa/Desktop/DECADE/Y6DESDM
export SWARP_DIR=/home/dhayaa/Desktop/DECADE/Y6DESDM/swarp-2.40.1/
export SRCEXT_DIR=/home/dhayaa/Desktop/DECADE/Y6DESDM/sextractor-2.24.4
export LSSTSTARSIM_DIR=/scratch/midway2/dhayaa/LSSTStars/


#DES-easyaccess
export DESREMOTE_RSYNC=desar2.cosmology.illinois.edu::ALLDESFiles/desarchive
export DESREMOTE_RSYNC_USER=dhayaa
export DES_RSYNC_PASSFILE=${HOME}/.desrsyncpass

export DECADEREMOTE_WGET=https://decade.ncsa.illinois.edu/deca_archive/

export EMAIL_ID=dhayaa@uchicago.edu
