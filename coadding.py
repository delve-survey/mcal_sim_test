import numpy as np
import yaml

from constants import MEDSCONF, PIFF_RUN
from files import get_band_info_file, make_dirs_for_file, get_swarp_flist_path
from des_info import add_extra_des_coadd_tile_info
import fitsio

TMP_DIR = os.environ['TMPDIR']

class MakeSwarpCoadds(object):
    """
    Class to create swarp coadds from simulated single epoch exposures.
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
                
                
        self.swarp_path = get_swarp_flist_path(output_meds_dir, MEDSCONF)
    
    def go(self):
        
        self._make_nwgint_files()
        self._make_filelists()
        self._make_coadds()
        self._make_detection_coadd()
            
    def _make_nwgint_files(self):
        '''
        Make null weight files using pixcorrect.
        Store them in temporary directory not project dir. 
        Mimics DESDM since they don't store these either.
        '''
        
        for band in bands:
            
            
            print("Creating nwgint images for %s band" % band)
            
            for src in self.info[band]['src_info']:
                
                args = {}
                
                pix_command = "$PIXCORRECT_DIR/bin/coadd_nwgint \
                                    -i red/D00233601_g_c50_r3650p01_immasked.fits.fz \
                                    -o nwgint/DES0130-4623_r5137p01_D00233601_g_c50_nwgint.fits \
                                    --headfile aux/DES0130-4623_r5137p01_D00233601_g_c50_scamp.ohead \
                                    --max_cols 50  \
                                    -v  \
                                    --interp_mask TRAIL,BPM  \
                                    --invalid_mask EDGE \
                                    --null_mask BPM,BADAMP,EDGEBLEED,EDGE,CRAY,SSXTALK,STREAK,TRAIL  \
                                    --block_size 5 \
                                    --tilename DES0130-4623  \
                                    --tileid 119590 \
                                    --me_wgt_keepmask STAR  \
                                    --hdupcfg $DESDM_CONFIG/Y6A1_v1_coadd_nwgint.config  \
                                    --streak_file $DESDM_CONFIG/Y3A2_v5_streaks_update-Y1234_FINALCUT_v1.fits" % args
                
                os.system(pix_command)
               
        
        return 1
        
    def _make_filelists(self):
        '''
        Takes list of exposures from yaml file and creates a text file
        where each line is filename/location. Will pass this into swarp.
        '''
        
        
        for band in bands:
            
            
            print("Writing list of files for %s band" % band)
            
            #We run swarp twice so make file list for both runs.
            #The sci and flx list are the same both times.
            for coadd_type in ['wgt', 'msk']:
                

                fname = os.path.join(self.swarp_path, self.tilename + '_swarp-band-%s-sci.list' % coadd_type)
                with open(fname, 'w') as f:
                    for src in self.info[band]['src_info']:
                        fobj.write("%s[0]\n" % se_info['image_path'])
                    
                print("Finished writing image list to %s" % fname)
            
            
                fname = os.path.join(self.swarp_path, self.tilename + '_swarp-band-%s-flx.list' % coadd_type)
                with open(fname, 'w') as f:
                    for src in self.info[band]['src_info']:
                        flux = 10. ** (0.4 * (30.0 - se_info['magzp'])) #Convert mag zeropoint to flux zp
                        fobj.write("%.16g\n" % flux)
                    
                    
            
            #However, the wgt files are different.
            #So running commands separately, outside the loop
            
            fname = os.path.join(self.swarp_path, self.tilename + '_swarp-band-wgt-wgt.list')
            with open(fname, 'w') as f:
                for src in self.info[band]['src_info']:
                    fobj.write("%s[2]\n" % se_info['weight_path'])
                    
            print("Finished writing flx scales (zeropoints) list to %s" % fname)
            
            
            fname = os.path.join(self.swarp_path, self.tilename + '_swarp-band-msk-wgt.list')
            with open(fname, 'w') as f:
                for src in self.info[band]['src_info']:
                    fobj.write("%s[1]\n" % se_info['weight_path'])
                    
            print("Finished writing flx scales (zeropoints) list to %s" % fname)
        
        return 1
        
    def _make_coadds(self):
        '''
        Runs swarp using all SE sources for a band
        '''
        
        args = {
                "list_prefix" : self.swarp_path + "/" + self.tilename,
                "out_prefix" : None,
                "RA": None,
                "DEC": None,
                }
        
        "$SWARP_DIR/swarp \
        @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-sci.list \
        -c config/Y6A1_v1_swarp.config \
        -WEIGHTOUT_NAME coadd/DES0130-4623_r5137p01_i_wgt.fits \
        -CENTER 22.632611,-46.386111 \
        -PIXEL_SCALE 0.263 \
        -FSCALE_DEFAULT @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-flx.list \
        -IMAGE_SIZE 10000,10000 \
        -IMAGEOUT_NAME coadd/DES0130-4623_r5137p01_i_sci.fits \
        -COMBINE_TYPE WEIGHTED \
        -WEIGHT_IMAGE @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-wgt.list \
        -NTHREADS 8 \
        -BLANK_BADPIXELS Y"
        
        for band in bands:
            
            coadd_file = self.info['image_path'].replace(TMP_DIR, self.output_meds_dir)
            make_dirs_for_file(coadd_file)
            
            #Is of the format "$DIR/{tilename}_{band}" without the .fits extension
            args['out_prefix'] = os.path.join(os.path.dirname(coadd_file), os.path.splitext(self.info['filename'])[0])
            
            
            #Get header from original coadd to get center in RA and DEC
            header = fitsio.read_header(self.info['image_path'], ext = 1)
            
            args['RA']  = header['CRVAL1']
            args['DEC'] = header['CRVAL2']
            
            
            swarp_command_wgt = "$SWARP_DIR/swarp \
                                        @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-sci.list \
                                        -c config/Y6A1_v1_swarp.config \
                                        -WEIGHTOUT_NAME coadd/DES0130-4623_r5137p01_i_wgt.fits \
                                        -CENTER 22.632611,-46.386111 \
                                        -PIXEL_SCALE 0.263 \
                                        -FSCALE_DEFAULT @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-flx.list \
                                        -IMAGE_SIZE 10000,10000 \
                                        -IMAGEOUT_NAME coadd/DES0130-4623_r5137p01_i_sci.fits \
                                        -COMBINE_TYPE WEIGHTED \
                                        -WEIGHT_IMAGE @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-wgt.list \
                                        -NTHREADS 8 \
                                        -BLANK_BADPIXELS Y" % args
            
            swarp_command_msk = "$SWARP_DIR/swarp \
                                        @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-msk-sci.list \
                                        -c config/Y6A1_v1_swarp.config \
                                        -WEIGHTOUT_NAME coadd/DES0130-4623_r5137p01_i_msk.fits \
                                        -CENTER 22.632611,-46.386111 \
                                        -PIXEL_SCALE 0.263 \
                                        -FSCALE_DEFAULT @list/band-swarp-msk/DES0130-4623_r5137p01_i_swarp-band-msk-flx.list \
                                        -IMAGE_SIZE 10000,10000 \
                                        -IMAGEOUT_NAME coadd/DES0130-4623_r5137p01_i_tmp-sci.fits \
                                        -COMBINE_TYPE WEIGHTED \
                                        -WEIGHT_IMAGE @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-msk-wgt.list \
                                        -NTHREADS 8 \
                                        -BLANK_BADPIXELS Y" % args
            
            
#             swarp_command_wgt = "$SWARP_DIR/swarp \
#                                         @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-sci.list \
#                                         -c config/Y6A1_v1_swarp.config \
#                                         -WEIGHTOUT_NAME coadd/DES0130-4623_r5137p01_i_wgt.fits \
#                                         -CENTER 22.632611,-46.386111 \
#                                         -PIXEL_SCALE 0.263 \
#                                         -FSCALE_DEFAULT @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-flx.list \
#                                         -IMAGE_SIZE 10000,10000 \
#                                         -IMAGEOUT_NAME coadd/DES0130-4623_r5137p01_i_sci.fits \
#                                         -COMBINE_TYPE WEIGHTED \
#                                         -WEIGHT_IMAGE @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-wgt-wgt.list \
#                                         -NTHREADS 8 \
#                                         -BLANK_BADPIXELS Y" % args
            
#             swarp_command_msk = "$SWARP_DIR/swarp \
#                                         @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-msk-sci.list \
#                                         -c config/Y6A1_v1_swarp.config \
#                                         -WEIGHTOUT_NAME coadd/DES0130-4623_r5137p01_i_msk.fits \
#                                         -CENTER 22.632611,-46.386111 \
#                                         -PIXEL_SCALE 0.263 \
#                                         -FSCALE_DEFAULT @list/band-swarp-msk/DES0130-4623_r5137p01_i_swarp-band-msk-flx.list \
#                                         -IMAGE_SIZE 10000,10000 \
#                                         -IMAGEOUT_NAME coadd/DES0130-4623_r5137p01_i_tmp-sci.fits \
#                                         -COMBINE_TYPE WEIGHTED \
#                                         -WEIGHT_IMAGE @list/band-swarp-wgt/DES0130-4623_r5137p01_i_swarp-band-msk-wgt.list \
#                                         -NTHREADS 8 \
#                                         -BLANK_BADPIXELS Y" % args
            
            
            os.system(swarp_command_wgt)
            print("Finished swarp wgt coadd for %s band" % band)
            
            os.system(swarp_command_msk)
            print("Finished swarp msk coadd for %s band" % band)
            
            
            args = {}
            command_assemble = "$DESPYFITS_DIR/coadd_assemble \
                                        --sci_file sci.fits \
                                        --wgt_file wgt.fits \
                                        --msk_file msk.fits  \
                                        --outname  coadd.fits \
                                        --xblock 10 \
                                        --yblock 3 \
                                        --maxcols 100 \
                                        --mincols 1 \
                                        --no-keep_sci_zeros \
                                        --magzero 30 \
                                        --tilename DES0130-4623 \
                                        --tileid 119590 \
                                        --interp_image MSK \
                                        --ydilate 3 " % args
            
            os.system(command_assemble)
            print("Finished assembling coadd for %s band" % band)
                        
        return 1
        
    def _make_detection_coadd():
        '''
        Combine swarp coadds to make detection coadd for
        Source Extractor to run on.
        '''
        
        args = {}
        swarp_command_wgt = "$SWARP_DIR/swarp sci[r].fits,sci[i].fits,sci[z].fits  \
                                    -c config/Y6A1_v1_swarp.config  \
                                    -WEIGHTOUT_NAME det_wgt.fits  \
                                    -CENTER 22.632611,-46.386111 \
                                    -RESAMPLE Y \
                                    -RESAMPLING_TYPE NEAREST  \
                                    -COPY_KEYWORDS BUNIT,TILENAME,TILEID \
                                    -PIXEL_SCALE 0.263 \
                                    -IMAGE_SIZE 10000,10000  \
                                    -IMAGEOUT_NAME det_sci.fits  \
                                    -COMBINE_TYPE AVERAGE \
                                    -WEIGHT_IMAGE wgt[r].fits,wgt[i].fits,wgt[z].fits  \
                                    -NTHREADS 8  -BLANK_BADPIXELS Y" % args
        
        swarp_command_msk = "$SWARP_DIR/swarp sci[r].fits,sci[i].fits,sci[z].fits  \
                                    -c config/Y6A1_v1_swarp.config  \
                                    -WEIGHTOUT_NAME det_msk.fits  \
                                    -CENTER 22.632611,-46.386111 \
                                    -RESAMPLE Y \
                                    -RESAMPLING_TYPE NEAREST  \
                                    -COPY_KEYWORDS BUNIT,TILENAME,TILEID \
                                    -PIXEL_SCALE 0.263 \
                                    -IMAGE_SIZE 10000,10000  \
                                    -IMAGEOUT_NAME det_tmpsci.fits  \
                                    -COMBINE_TYPE AVERAGE \
                                    -WEIGHT_IMAGE msk[r].fits,msk[i].fits,msk[z].fits  \
                                    -NTHREADS 8  -BLANK_BADPIXELS Y" % args
        
        os.system(swarp_command_wgt)
        print("Finished swarp wgt coadd for %s band" % band)

        os.system(swarp_command_msk)
        print("Finished swarp msk coadd for %s band" % band)
        
        args = {}
        command_assemble = "$DESPYFITS_DIR/coadd_assemble \
                                        --sci_file det_sci.fits  \
                                        --wgt_file det_wgt.fits \
                                        --msk_file det_msk.fits  \
                                        --band det  \
                                        --outname NAME_det.fits  \
                                        --xblock 10 \
                                        --yblock 3 \
                                        --maxcols 100 \
                                        --mincols 1 \
                                        --no-keep_sci_zeros \
                                        --magzero 30 \
                                        --tilename DES0130-4623 \
                                        --tileid 119590 \
                                        --interp_image MSK \
                                        --ydilate 3" % args
            
        os.system(command_assemble)
        print("Finished assembling coadd")
        
        return 1
