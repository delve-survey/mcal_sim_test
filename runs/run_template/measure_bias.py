from astropy.io import fits
import numpy as np
import os
from tqdm import tqdm
import sys

from render_jobs import tilenames

print('--------------------------')
print('USING TILES:')
print('--------------------------')
for t in tilenames:
    print(t)
print('--------------------------')

N = len(tilenames)

name     = os.path.basename(os.path.dirname(__file__))
MCAL_DIR = os.environ['MCAL_DIR']
PATH     = MCAL_DIR + '/' + name
print('GETTING MCAL FILES FROM:')
print(PATH)

tiles_list = [[x] for x in tilenames]
tiles_list.append(tilenames)

#Loop over tiles (and also a combination of all tiles together)
#and get m and c for objects in each tile
for tiles in tiles_list:
    N = len(tiles)
    
    
    print("--------------------------------------------------------")
    print("TILES:", tiles)
    print("--------------------------------------------------------")

    gplus  = [fits.open(PATH + '/metacal_%s_seed0_gplus.fits'%t) for t in tiles]
    gminus = [fits.open(PATH + '/metacal_%s_seed0_gminus.fits'%t) for t in tiles]

    shear = '_noshear'
    Flags_plus = np.concatenate([gplus[i][1].data['mcal_flags'] for i in range(N)])
    SNR_plus   = np.concatenate([gplus[i][1].data['mcal_s2n' + shear] for i in range(N)])
    T_plus     = np.concatenate([gplus[i][1].data['mcal_T_ratio' + shear] for i in range(N)])
    
    Mask_plus  = (Flags_plus == 0) & (SNR_plus > 10) & (T_plus > 0.5)
    
    print("FRACTION KEPT (PLUS) :", Mask_plus.sum()/Mask_plus.size)
    
    Flags_minus = np.concatenate([gminus[i][1].data['mcal_flags'] for i in range(N)])
    SNR_minus   = np.concatenate([gminus[i][1].data['mcal_s2n' + shear] for i in range(N)])
    T_minus     = np.concatenate([gminus[i][1].data['mcal_T_ratio' + shear] for i in range(N)])
    
    Mask_minus  = (Flags_minus == 0) & (SNR_minus > 10) & (T_minus > 0.5)
    print("FRACTION KEPT (MINUS):",Mask_minus.sum()/Mask_minus.size)
    
    e1_plus_array = np.concatenate([gplus[i][1].data['mcal_g_noshear'][:,0] for i in range(N)])[Mask_plus]
    e2_plus_array = np.concatenate([gplus[i][1].data['mcal_g_noshear'][:,1] for i in range(N)])[Mask_plus]

    e1_minus_array = np.concatenate([gminus[i][1].data['mcal_g_noshear'][:,0] for i in range(N)])[Mask_minus]
    e2_minus_array = np.concatenate([gminus[i][1].data['mcal_g_noshear'][:,1] for i in range(N)])[Mask_minus]


    R11_plus_array = np.concatenate([(gplus[i][1].data['mcal_g_1p'][:,0] - gplus[i][1].data['mcal_g_1m'][:,0])/0.02 for i in range(N)])[Mask_plus]
    R22_plus_array = np.concatenate([(gplus[i][1].data['mcal_g_2p'][:,1] - gplus[i][1].data['mcal_g_2m'][:,1])/0.02 for i in range(N)])[Mask_plus]


    R11_minus_array = np.concatenate([(gminus[i][1].data['mcal_g_1p'][:,0] - gminus[i][1].data['mcal_g_1m'][:,0])/0.02 for i in range(N)])[Mask_minus]
    R22_minus_array = np.concatenate([(gminus[i][1].data['mcal_g_2p'][:,1] - gminus[i][1].data['mcal_g_2m'][:,1])/0.02 for i in range(N)])[Mask_minus]

    del gplus, gminus

    nBoot = 1000

    e1_plus  = np.zeros(nBoot)
    e1_minus = np.zeros(nBoot)
    e2_plus  = np.zeros(nBoot)
    e2_minus = np.zeros(nBoot)

    R11_plus  = np.zeros(nBoot)
    R11_minus = np.zeros(nBoot)
    R22_plus  = np.zeros(nBoot)
    R22_minus = np.zeros(nBoot)

    
    print("SIZE : ", e1_plus_array.shape)


    for i in tqdm(range(nBoot)):

        index_plus  = np.random.randint(0, len(e1_plus_array),  len(e1_plus_array))
        index_minus = np.random.randint(0, len(e1_minus_array), len(e1_minus_array))

        e1_plus[i]  = np.mean(e1_plus_array[index_plus])
        e1_minus[i] = np.mean(e1_minus_array[index_minus])
        e2_plus[i]  = np.mean(e2_plus_array[index_plus])
        e2_minus[i] = np.mean(e2_minus_array[index_minus])

        R11_plus[i]   = np.mean(R11_plus_array[index_plus])
        R11_minus[i]  = np.mean(R11_minus_array[index_minus])
        R22_plus[i]   = np.mean(R22_plus_array[index_plus])
        R22_minus[i]  = np.mean(R22_minus_array[index_minus])


    R11 = (R11_plus + R11_minus)/2
    R22 = (R22_plus + R22_minus)/2

    g1 = (e1_plus - e1_minus)/(2*R11)
    g2 = (e2_plus + e2_minus)/(2*R22)

    m = g1/0.02 - 1
    c = g2

    print('\n---------------------------------\n')
    print("----------- COMBINED ---------------")
    print('size      N = %d'%e1_plus_array.shape)
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)))
    
    print('\n---------------------------------\n')
    
    R11 = R11_plus
    R22 = R22_plus

    g1 = e1_plus/R11
    g2 = e2_plus/R22
    
    m = g1/0.02 - 1
    c = g2

    print("----------- PLUS ---------------")
    print('size      N = %d'%e1_plus_array.shape)
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]\n\n'%(np.mean(c/1e-5), 3*np.std(c/1e-5)))
    
    R11 = R11_minus
    R22 = R22_minus

    g1 = e1_minus/R11
    g2 = e2_minus/R22
    
    m = g1/(-0.02) - 1
    c = g2

    print("----------- MINUS ---------------")
    print('size      N = %d'%e1_minus_array.shape)
    print('recovered m = %1.3f +/- %1.3f [1e-3, 3-sigma]'%(np.mean(m/1e-3), 3*np.std(m/1e-3)))
    print('recovered c = %1.3f +/- %1.3f [1e-5, 3-sigma]'%(np.mean(c/1e-5), 3*np.std(c/1e-5)))
          
          
    print("------------------------------------------------------------")
