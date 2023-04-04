'''
Routines for moving files around, deleting remaining files
and saving relevant quantities.
'''

import shutil
import yaml
import os

def finalize_files(tilename, bands, output_desdata):
    
    for b in bands:
        move_SrcExtractor_cat(tilename, b, output_desdata)
        move_meds(tilename, b, output_desdata)
        
    move_metacal_cat(tilename, output_desdata)
    #clean_up_directory()
        

#Helper functions to run the above cleanup/re-organization code
def move_SrcExtractor_cat(tile, band, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'tile' : tile,
            'band' : band}
    
    #This is always going to be path in our runs so just sort of hardcode this assumption
    config_path = output_desdata + '/simple_des_y3_sims/y3v02/band_info_files/%(tile)s_%(band)s_info.yaml'%args
    with open(config_path, 'r') as fp:
        band_info = yaml.load(fp, Loader=yaml.Loader)

    cat_path = band_info['cat_path'].replace(os.environ['TMPDIR'], output_desdata)

    new_path = os.environ['MCAL_DIR'] + "/%(name)s/SrcExtractor_%(tile)s_g%(mode)s_%(band)s-cat.fits" % args

    print(cat_path, new_path)
    shutil.move(cat_path, new_path)
        
    return True


def move_metacal_cat(tile, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'tile' : tile}
    
    cat_path = output_desdata + "/metacal/y3v02/%(tile)s_metacal.fits" % args 
    new_path = os.environ['MCAL_DIR'] + "/%(name)s/metacal_%(tile)s_g%(mode)s.fits" % args

    print(cat_path, new_path)
    shutil.move(cat_path, new_path)
         
    return True

def move_meds(tile, band, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'tile' : tile,
            'band' : band}
    
    meds_path = output_desdata + "/meds/y3v02/%(tile)s/%(tile)s_%(band)s_meds-y3v02.fits.fz" % args 
    new_path  = os.environ['MCAL_DIR'] + "/%(name)s/meds_%(tile)s_g%(mode)s_%(band)s-y3v02.fits.fz" % args

    print(meds_path, new_path)
    shutil.move(meds_path, new_path)
        
    return True

def clean_up_directory():
    shutil.rmtree("$output")
    
    return True
