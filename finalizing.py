'''
Routines for moving files around, deleting remaining files
and saving relevant quantities.
'''

import shutil
import yaml
import os

def finalize_files(tile, bands, seed):
    
    for b in bands:
        move_SrcExtractor_cat(tile, b, seed)
        move_meds(tile, b, seed)
        
    move_metacal_cat(tile, seed)
    clean_up_directory()
        

#Helper functions to run the above cleanup/re-organization code
def move_SrcExtractor_cat(tile, band, seed):
    
    for mode in ['plus', 'minus']:
        args = {'name' : name,
                'seed' : seed,
                'mode' : mode,
                'tile' : tile,
                'band' : band}

        #This is always going to be path in our runs so just sort of hardcode this assumption
        config_path = os.environ['PREP_DIR'] + '/%(name)s/outputs_%(tile)s_seed%(seed)d_g%(mode)s/simple_des_y3_sims/y3v02/band_info_files/%(tile)s_%(band)s_info.yaml'%args
        with open(config_path, 'r') as fp:
            band_info = yaml.load(fp, Loader=yaml.Loader)

        cat_path = band_info['cat_path'].replace(os.environ['TMPDIR'], 
                                                 os.environ['PREP_DIR'] + '/%(name)s/outputs_%(tile)s_seed%(seed)d_g%(mode)s/'%args)
        
        new_path = "$MCAL_DIR/%(name)s/SrcExtractor_%(tile)s_seed%(seed)d_g%(mode)s_%(band)-cat.fits" % args
        
        shutil.move(cat_path, new_path)
        
    return True


def move_metacal_cat(tile, seed):
    
    for mode in ['plus', 'minus']:
        args = {'name' : name,
                'seed' : seed,
                'mode' : mode,
                'tile' : tile,
                'band' : band}
    
        cat_path = "$output/metacal/y3v02/%(tile)s_metacal.fits" 
        new_path = "$MCAL_DIR/%(name)s/metacal_%s(tile)_seed%(seed)d_g%(mode).fits"
    
    return True

def move_meds(tile, band, seed):
    
    for mode in ['plus', 'minus']:
        args = {'name' : name,
                'seed' : seed,
                'mode' : mode,
                'tile' : tile,
                'band' : band}
    
        meds_path = "$output/meds/y3v02/%(tile)s/%(tile)s_%(band)s_meds-y3v02.fits.fz" 
        new_path  = "$MCAL_DIR/%(name)s/meds_%(tile)s_g%(mode)s_%(band)s-y3v02.fits.fz"
        
        shutil.move(meds_path, new_path)
        
    return True

def clean_up_directory():
    shutil.rmtree("$output")
    
    return True