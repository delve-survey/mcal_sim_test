import numpy as np
import yaml
import os
from constants import MEDSCONF, PIFF_RUN
from files import (get_band_info_file, make_dirs_for_file, 
                   get_swarp_files_path, get_nwgint_path)
from des_info import add_extra_des_coadd_tile_info
import fitsio
import shutil
from constants import MEDSCONF, MAGZP_REF

TMP_DIR = os.environ['TMPDIR']

class MakeSrcExtractorCat(object):
    
    """
    Class to get sourceExtractor catalog from detection coadd
    """
    
    def __init__(self, *, tilename, bands, output_meds_dir, config, n_files=None):
        
        self.output_meds_dir = output_meds_dir
        self.tilename = tilename
        self.bands = bands
        
        self.info = {}
        for band in bands:
            fname = get_band_info_file(
                meds_dir=self.output_meds_dir,
                medsconf=MEDSCONF,
                tilename=self.tilename,
                band=band)
            with open(fname, 'r') as fp:
                self.info[band] = yaml.load(fp, Loader=yaml.Loader)
                
        self.swarp_path = get_swarp_files_path(meds_dir=output_meds_dir, medsconf=MEDSCONF)
        
    def run(self):
        
        for band in self.bands:
            coadd_file = self.info[band]['image_path'].replace(TMP_DIR, self.output_meds_dir).replace('.fz', '')
            seg_file   = self.info[band]['seg_path'].replace(TMP_DIR, self.output_meds_dir)
            cat_file   = self.info[band]['cat_path'].replace(TMP_DIR, self.output_meds_dir)
            
            make_dirs_for_file(cat_file)
            make_dirs_for_file(seg_file)
            
            #Copy segmentation file over
            #Needed for meds-making code
            shutil.copy(self.info[band]['seg_path'], seg_file)

            #Setup args
            args = {'DET_COADD' : os.path.join(self.swarp_path, self.tilename + '_det.fits'),
                    'TILENAME'  : self.tilename,
                    'BAND'      : band,
                    'COADD'     : coadd_file,
                    'SEG'       : seg_file,
                    'CAT'       : cat_file}


            srcextractor_command = "$SRCEXT_DIR/src/sex \
                                        %(DET_COADD)s[0],%(COADD)s[0] \
                                        -c $DESDM_CONFIG/Y6A1_v1_sex.config \
                                        -CHECKIMAGE_TYPE SEGMENTATION \
                                        -CHECKIMAGE_NAME %(SEG)s \
                                        -PARAMETERS_NAME $DESDM_CONFIG/Y6A1_v1_srcex.param_diskonly_old \
                                        -MAG_ZEROPOINT 30 \
                                        -FILTER_NAME $DESDM_CONFIG/Y6A1_v1_gauss_3.0_7x7.conv \
                                        -CATALOG_NAME %(CAT)s \
                                        -FLAG_IMAGE %(COADD)s[1] \
                                        -STARNNW_NAME $DESDM_CONFIG/Y6A1_v1_srcex.nnw \
                                        -WEIGHT_IMAGE %(DET_COADD)s[2],%(COADD)s[2] \
                                        -DEBLEND_MINCONT 0.001 \
                                        -DEBLEND_NTHRESH 64 \
                                        -DETECT_THRESH 0.8 \
                                        -ANALYSIS_THRESH 0.8" % args
            
            os.system(srcextractor_command)
            print("Finished srcextractor run for %s band" % band)
            