'''
Routines for moving files around, deleting remaining files
and saving relevant quantities.
'''

import shutil
import yaml
import os

def initialize_files(tilename, bands, output_desdata, config):
    
    for b in bands:
        rewrite_band_infoing(tilename, b, output_desdata)
        
    move_prepdir_files(tilename)
    move_tmpdir_files(tilename)        

#Helper functions to run the above cleanup/re-organization code
def move_tmpdir_files(tile, band, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'tile' : tile,
            'band' : band}
    
    old_path = os.environ['SHEAR_TMPDIR'] + "/%(name)s/" % args
    new_path = os.environ['TMPDIR'] + "/%(name)s/" % args

    print(old_path, new_path)
    shutil.move(old_path, new_path)
        
    return True

#Helper functions to run the above cleanup/re-organization code
def move_prepdir_files(tile, band, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'basename': os.path.basename(output_desdata),
            'tile' : tile,
            'band' : band}
    
    old_path = os.environ['PREPDIR'] + "/%(name)s/%(basename)s" % args
    new_path = os.environ['TMPDIR'] + "/%(name)s/%(basename)s" % args

    print(old_path, new_path)
    shutil.move(old_path, new_path)
        
    return True


def rewrite_band_infoing(tile, band, output_desdata):
    
    def rewrite(X, a, b):
        
        for k in X.keys():
            
            if isinstance(X[k], str):
                X[k] = X[k].replace(a, b)
         
            if isinstance(X[k], dict):
                X[k] = rewrite(X, a, b)
                
        return X
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'mode' : 'plus' if 'plus' in os.path.basename(output_desdata) else 'minus',
            'tile' : tile,
            'band' : band}
    
    #This is always going to be path in our runs so just sort of hardcode this assumption
    band_path = output_desdata + '/simple_des_y3_sims/y3v02/band_info_files/%(tile)s_%(band)s_info.yaml'%args
    with open(band_path, 'r') as fp:
        band_info = yaml.load(fp, Loader=yaml.Loader)
    
    old_path = os.environ['SHEAR_TMPDIR'] + "/%(name)s/" % args
    new_path = os.environ['TMPDIR'] + "/%(name)s/" % args
    
    band_info_new = rewrite(band_info, old_path, new_path)
    
    with open(band_path, 'w') as fp:
        yaml.dump(band_info, fp)
         
    return True