# mcal_sim_test
Testing metacal measurements on simulated images

## Basic instructions:
1. create a conda environment using ```conda env create --file environment.yml```

2. activate the environment by doing ```conda activate shear```

3. that's probably not a complete set of all packages you need; in particular you should download desmeds from https://github.com/esheldon/desmeds (because the conda version doesn't have all we need); you might or might not need to install desmeds manually following the instructions in that link; you'll also need easyaccess and other simpler conda-installable packages

4. to run the basic prep stage, do for instance: ```python run_sims.py prep --tilename="DES0544-2249" --bands="riz" --output-desdata="outputs-DES0544-2249"```
