import numpy as np
import yaml

from constants import MEDSCONF, PIFF_RUN
from files import (get_band_info_file, make_dirs_for_file, 
                   get_swarp_files_path, get_nwgint_path)
from des_info import add_extra_des_coadd_tile_info
import fitsio
from constants import MEDSCONF, MAGZP_REF

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
                
                
        self.swarp_path = get_swarp_files_path(output_meds_dir, MEDSCONF)
        
        self.nwgint_flists = {b:[] for b in bands}
    
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
        
        for band in self.bands:
            
            #Get header from original coadd to get tileid (but also tilename since we can)
            header = fitsio.read_header(self.info[band]['image_path'], ext = 1)
           
            args = {}
            
            args['TILENAME'] = header['TILENAME']
            args['TILEID']   = header['TILEID']
            
            #Get base output and make dirs if needed
            out_basepath = get_nwgint_path(self.output_meds_dir, MEDSCONF, band)
            make_dirs_for_file(out_basepath)
            
            print("Creating nwgint images for %s band" % band)
            
            for src in self.info[band]['src_info']:
                
                args['IMAGE_PATH'] = src['image_path']
                args['HEAD_PATH']  = src['head_path']
                args['OUT_PATH']   = os.join(out_basepath, src['filename'].replace('immasked', 'nwgint'))
                                
                #add nwgint info back into THIS specific info dict
                src['nwgint_path'] = args['OUT_PATH']
                
                pix_command = "$PIXCORRECT_DIR/bin/coadd_nwgint \
                                    -i %(IMAGE_PATH)s \
                                    -o %(OUT_PATH)s \
                                    --headfile %(HEAD_PATH)s \
                                    --max_cols 50  \
                                    -v  \
                                    --interp_mask TRAIL,BPM  \
                                    --invalid_mask EDGE \
                                    --null_mask BPM,BADAMP,EDGEBLEED,EDGE,CRAY,SSXTALK,STREAK,TRAIL  \
                                    --block_size 5 \
                                    --tilename %(TILENAME)s  \
                                    --tileid %(TILEID)s \
                                    --me_wgt_keepmask STAR  \
                                    --hdupcfg $DESDM_CONFIG/Y6A1_v1_coadd_nwgint.config  \
                                    --streak_file $DESDM_CONFIG/Y3A2_v6_streaks_update-Y1234_FINALCUT_v2.fits" % args
                
#                 pix_command = "$PIXCORRECT_DIR/bin/coadd_nwgint \
#                                     -i red/D00233601_g_c50_r3650p01_immasked.fits.fz \
#                                     -o nwgint/DES0130-4623_r5137p01_D00233601_g_c50_nwgint.fits \
#                                     --headfile aux/DES0130-4623_r5137p01_D00233601_g_c50_scamp.ohead \
#                                     --max_cols 50  \
#                                     -v  \
#                                     --interp_mask TRAIL,BPM  \
#                                     --invalid_mask EDGE \
#                                     --null_mask BPM,BADAMP,EDGEBLEED,EDGE,CRAY,SSXTALK,STREAK,TRAIL  \
#                                     --block_size 5 \
#                                     --tilename DES0130-4623  \
#                                     --tileid 119590 \
#                                     --me_wgt_keepmask STAR  \
#                                     --hdupcfg $DESDM_CONFIG/Y6A1_v1_coadd_nwgint.config  \
#                                     --streak_file $DESDM_CONFIG/Y3A2_v5_streaks_update-Y1234_FINALCUT_v1.fits" % args
                
                os.system(pix_command)
               
        
        return 1
        
    def _make_filelists(self):
        '''
        Takes list of exposures from yaml file and creates a text file
        where each line is filename/location. Will pass this into swarp.
        '''
        
        
        for band in self.bands:
            
            
            print("Writing list of files for %s band" % band)
            
            #We run swarp twice so make file list for both runs.
            #The sci and flx list are the same both times.
            for coadd_type in ['wgt', 'msk']:
                

                fname = os.path.join(self.swarp_path, self.tilename + '_swarp-band-%s-sci.list' % coadd_type)
                with open(fname, 'w') as f:
                    for src in self.info[band]['src_info']:
                        fobj.write("%s[0]\n" % src['nwgint_path'])
                    
                print("Finished writing image list to %s" % fname)
            
            
                fname = os.path.join(self.swarp_path, self.tilename + '_swarp-band-%s-flx.list' % coadd_type)
                with open(fname, 'w') as f:
                    for src in self.info[band]['src_info']:
                        flux = 10. ** (0.4 * (30.0 - src['magzp'])) #Convert mag zeropoint to flux zp
                        fobj.write("%.16g\n" % flux)
                    
                    
            
            #However, the wgt files are different.
            #So running commands separately, outside the loop
            
            fname = os.path.join(self.swarp_path, self.tilename + '_swarp-band-wgt-wgt.list')
            with open(fname, 'w') as f:
                for src in self.info[band]['src_info']:
                    fobj.write("%s[2]\n" % src['weight_path'])
                    
            print("Finished writing flx scales (zeropoints) list to %s" % fname)
            
            
            fname = os.path.join(self.swarp_path, self.tilename + '_swarp-band-msk-wgt.list')
            with open(fname, 'w') as f:
                for src in self.info[band]['src_info']:
                    fobj.write("%s[1]\n" % src['weight_path'])
                    
            print("Finished writing flx scales (zeropoints) list to %s" % fname)
        
        return 1
        
    def _make_coadds(self):
        '''
        Runs swarp using all SE sources for a band
        '''
        
        #Use first band to get the RA and DEC for the tile center.
        #center is same for all bands
        
        #Get header from original coadd to get center in RA and DEC
        header = fitsio.read_header(self.info[self.bands[0]]['image_path'], ext = 1)


        args = {
                "list_prefix" : self.swarp_path + "/" + self.tilename,
                "out_prefix"  : None,
                "RA" : header['CRVAL1'],
                "DEC": header['CRVAL2'],
                "TILENAME" : header['TILENAME'],
                "TILEID"   : header['TILEID']
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
        
        for band in self.bands:
            
            coadd_file = self.info[band]['image_path'].replace(TMP_DIR, self.output_meds_dir)
            make_dirs_for_file(coadd_file)
            
            #Is of the format "$DIR/{tilename}_{band}" without the .fits.fz extension
            args['out_prefix'] = coadd_file.replace('.fits.fz', '')
            
            
            swarp_command_wgt = "$SWARP_DIR/src/swarp \
                                        @%(list_prefix)s-wgt-sci.list \
                                        -c $DESDM_CONFIG/Y6A1_v1_swarp.config \
                                        -WEIGHTOUT_NAME %(out_prefix)s_wgt.fits \
                                        -CENTER %(RA)0.6f, %(DEC)0.6f \
                                        -PIXEL_SCALE 0.263 \
                                        -FSCALE_DEFAULT @%(list_prefix)s-wgt-flx.list \
                                        -IMAGE_SIZE 10000,10000 \
                                        -IMAGEOUT_NAME %(out_prefix)s_sci.fits \
                                        -COMBINE_TYPE WEIGHTED \
                                        -WEIGHT_IMAGE @%(list_prefix)s-wgt-wgt.list \
                                        -NTHREADS 8 \
                                        -BLANK_BADPIXELS Y" % args
            
            
            swarp_command_msk = "$SWARP_DIR/src/swarp \
                                        @%(list_prefix)s-msk-sci.list \
                                        -c $DESDM_CONFIG/Y6A1_v1_swarp.config \
                                        -WEIGHTOUT_NAME %(out_prefix)s_msk.fits \
                                        -CENTER %(RA)0.6f, %(DEC)0.6f \
                                        -PIXEL_SCALE 0.263 \
                                        -FSCALE_DEFAULT @%(list_prefix)s-msk-flx.list \
                                        -IMAGE_SIZE 10000,10000 \
                                        -IMAGEOUT_NAME %(out_prefix)s_tmp-sci.fits \
                                        -COMBINE_TYPE WEIGHTED \
                                        -WEIGHT_IMAGE @%(list_prefix)s-msk-wgt.list \
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
            command_assemble = "$DESPYFITS_DIR/bin/coadd_assemble \
                                        --sci_file %(out_prefix)s_sci.fits \
                                        --wgt_file %(out_prefix)s_wgt.fits \
                                        --msk_file %(out_prefix)s_msk.fits  \
                                        --outname  %(out_prefix)s.fits \
                                        --xblock 10 \
                                        --yblock 3 \
                                        --maxcols 100 \
                                        --mincols 1 \
                                        --no-keep_sci_zeros \
                                        --magzero 30 \
                                        --tilename %(TILENAME)s \
                                        --tileid %(TILEID)s \
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
        
        
        #Get header from original coadd to get center in RA and DEC
        header = fitsio.read_header(self.info[self.bands[0]]['image_path'], ext = 1)


        sci_paths = ['%s_sci.fits'%self.info[b]['image_path'].replace(TMP_DIR, self.output_meds_dir).replace('.fits.fz', '') for b in self.bands]
        wgt_paths = ['%s_wgt.fits'%self.info[b]['image_path'].replace(TMP_DIR, self.output_meds_dir).replace('.fits.fz', '') for b in self.bands]
        msk_paths = ['%s_msk.fits'%self.info[b]['image_path'].replace(TMP_DIR, self.output_meds_dir).replace('.fits.fz', '') for b in self.bands]
        
        
        args = {
                "sci_paths"   : ','.join(sci_paths), #Convert from list to single string with ',' separator
                "wgt_paths"   : ','.join(wgt_paths),
                "msk_paths"   : ','.join(msk_paths),
                "out_prefix"  : self.swarp_path,
                "RA" : header['CRVAL1'],
                "DEC": header['CRVAL2'],
                "TILENAME" : header['TILENAME'],
                "TILEID"   : header['TILEID']
                }
        
        swarp_command_wgt = "$SWARP_DIR/src/swarp s%(sci_paths)s  \
                                    -c $DESDM_CONFIG/Y6A1_v1_swarp.config  \
                                    -WEIGHTOUT_NAME %(out_prefix)s/%(TILENAME)s_det_wgt.fits  \
                                    -CENTER %(RA)0.6f, %(DEC)0.6f \
                                    -RESAMPLE Y \
                                    -RESAMPLING_TYPE NEAREST  \
                                    -COPY_KEYWORDS BUNIT,TILENAME,TILEID \
                                    -PIXEL_SCALE 0.263 \
                                    -IMAGE_SIZE 10000,10000  \
                                    -IMAGEOUT_NAME %(out_prefix)s/%(TILENAME)s_det_sci.fits  \
                                    -COMBINE_TYPE AVERAGE \
                                    -WEIGHT_IMAGE %(wgt_paths)s  \
                                    -NTHREADS 8  \
                                    -BLANK_BADPIXELS Y" % args
        
        swarp_command_msk = "$SWARP_DIR/src/swarp s%(sci_paths)s  \
                                    -c $DESDM_CONFIG/Y6A1_v1_swarp.config  \
                                    -WEIGHTOUT_NAME %(out_prefix)s/%(TILENAME)s_det_msk.fits  \
                                    -CENTER %(RA)0.6f, %(DEC)0.6f \
                                    -RESAMPLE Y \
                                    -RESAMPLING_TYPE NEAREST  \
                                    -COPY_KEYWORDS BUNIT,TILENAME,TILEID \
                                    -PIXEL_SCALE 0.263 \
                                    -IMAGE_SIZE 10000,10000  \
                                    -IMAGEOUT_NAME %(out_prefix)s/%(TILENAME)s_det_tmpsci.fits  \
                                    -COMBINE_TYPE AVERAGE \
                                    -WEIGHT_IMAGE %(msk_paths)s  \
                                    -NTHREADS 8  \
                                    -BLANK_BADPIXELS Y" % args
        
#         swarp_command_wgt = "$SWARP_DIR/src/swarp sci[r].fits,sci[i].fits,sci[z].fits  \
#                                     -c config/Y6A1_v1_swarp.config  \
#                                     -WEIGHTOUT_NAME det_wgt.fits  \
#                                     -CENTER 22.632611,-46.386111 \
#                                     -RESAMPLE Y \
#                                     -RESAMPLING_TYPE NEAREST  \
#                                     -COPY_KEYWORDS BUNIT,TILENAME,TILEID \
#                                     -PIXEL_SCALE 0.263 \
#                                     -IMAGE_SIZE 10000,10000  \
#                                     -IMAGEOUT_NAME det_sci.fits  \
#                                     -COMBINE_TYPE AVERAGE \
#                                     -WEIGHT_IMAGE wgt[r].fits,wgt[i].fits,wgt[z].fits  \
#                                     -NTHREADS 8  -BLANK_BADPIXELS Y" % args
        
#         swarp_command_msk = "$SWARP_DIR/src/swarp sci[r].fits,sci[i].fits,sci[z].fits  \
#                                     -c config/Y6A1_v1_swarp.config  \
#                                     -WEIGHTOUT_NAME det_msk.fits  \
#                                     -CENTER 22.632611,-46.386111 \
#                                     -RESAMPLE Y \
#                                     -RESAMPLING_TYPE NEAREST  \
#                                     -COPY_KEYWORDS BUNIT,TILENAME,TILEID \
#                                     -PIXEL_SCALE 0.263 \
#                                     -IMAGE_SIZE 10000,10000  \
#                                     -IMAGEOUT_NAME det_tmpsci.fits  \
#                                     -COMBINE_TYPE AVERAGE \
#                                     -WEIGHT_IMAGE msk[r].fits,msk[i].fits,msk[z].fits  \
#                                     -NTHREADS 8  -BLANK_BADPIXELS Y" % args
        
        os.system(swarp_command_wgt)
        print("Finished swarp wgt coadd for %s band" % band)

        os.system(swarp_command_msk)
        print("Finished swarp msk coadd for %s band" % band)
        
        args = {}
        command_assemble = "$DESPYFITS_DIR/bin/coadd_assemble \
                                        --sci_file %(out_prefix)s/%(TILENAME)s_det_sci.fits  \
                                        --wgt_file %(out_prefix)s/%(TILENAME)s_det_wgt.fits \
                                        --msk_file %(out_prefix)s/%(TILENAME)s_det_msk.fits  \
                                        --band det  \
                                        --outname %(out_prefix)s/%(TILENAME)s_det.fits  \
                                        --xblock 10 \
                                        --yblock 3 \
                                        --maxcols 100 \
                                        --mincols 1 \
                                        --no-keep_sci_zeros \
                                        --magzero 30 \
                                        --tilename %(TILENAME)s \
                                        --tileid %(TILEID)s \
                                        --interp_image MSK \
                                        --ydilate 3" % args
            
        os.system(command_assemble)
        print("Finished assembling coadd")
        
#         os.system(r'rm %(out_prefix)s/%(TILENAME)s_det_tmpsci.fits' % args)
#         print("Finished removing temp sci file)
        
        return 1