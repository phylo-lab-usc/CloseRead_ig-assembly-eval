source ${7}
conda config --append envs_dirs ${8}
conda activate ${8}/IGdetective
conda env list
conda info --envs

export PATH=${8}/IGdetective/bin/python:$PATH
export PATH=$(echo $PATH | sed -e 's|/spack/conda/miniforge3/24.3.0/bin:||g')

echo "PATH: $PATH"
which python

${8}/IGdetective/bin/python ${6}/run_iterative_igdetective.py $1 $2

touch ${5}/igGene/${3}.${4}.txt
