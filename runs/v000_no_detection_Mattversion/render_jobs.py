import os
import jinja2

tilenames = [
    'DES0544-2249',
    'DES2122+0001',
    'DES2122-0041',
    'DES2122+0043',
    'DES2122-0124',
    'DES2122+0126',
    'DES2122-0207',
    'DES2122+0209',
    'DES2125+0001',
    'DES2125-0041',
    'DES2125+0043',
    'DES2125-0124',
    'DES2125+0126',
    'DES2125-0207',
    'DES2125+0209']
# 'DES2128+0001',
# 'DES2128-0041',
# 'DES2128+0043',
# 'DES2128-0124',
# 'DES2128+0126',
# 'DES2128-0207',
# 'DES2128+0209',
# 'DES2131+0001',
# 'DES2131-0041',
# 'DES2131+0043',
# 'DES2131-0124',
# 'DES2131+0126',
# 'DES2131-0207',
# 'DES2131+0209',
# 'DES2134+0001',
# 'DES2134-0041',
# 'DES2134+0043',
# 'DES2134-0124',
# 'DES2134+0126',
# 'DES2134-0207',
# 'DES2134+0209',
# 'DES2137+0001',
# 'DES2137-0041',
# 'DES2137+0043',
# 'DES2137-0124',
# 'DES2137+0126',
# 'DES2137+0209']

with open('job.sh.temp', 'r') as fp:
    tmp = jinja2.Template(fp.read())

for i, tilename in enumerate(tilenames):
    with open('job_%s.sh' % tilename, 'w') as fp:
        fp.write(tmp.render(
            tilename=tilename, seed_galsim=i+1, seed_mcal=i+1000))
    os.system('chmod u+x job_%s.sh' % tilename)