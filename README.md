
<!-- README.md is generated from README.Rmd. Please edit that file -->

# test-opti-random

<!-- badges: start -->
<!-- badges: end -->

Testing the changes Sarah Westcott made to OptiFit in the mothur 1.47.0
branch.

## Setup

On GreatLakes, clone the repos and compile mothur version 1.47.0

``` bash
git clone https://github.com/kelly-sovacool/test-opti-random
git clone https://github.com/mothur/mothur
cd mothur
git checkout 1.47.0
module load gcc/4.8.5   
module load boost/1.70.0
make -j 4 -f Makefile_cluster USEBOOST=yes BOOST_LIBRARY_DIR=${BOOST_LIB} BOOST_INCLUDE_DIR=${BOOST_INCLUDE} USEHDF5=no install
cd ..
cp mothur/mothur test-opti-random/bin/mothur-1.47.0/
```

Install other dependencies and activate the conda environment:

``` bash
cd test-opti-random
conda env create -f config/env.simple.yml
conda activate optifit
```

Run the workflow:

    locally:
    ```bash
    snakemake -j 4
    ```

    or on the cluster: (modify your email and account name first)
    ```
    sbatch code/submit_KLS.sh
    ```

## Results
