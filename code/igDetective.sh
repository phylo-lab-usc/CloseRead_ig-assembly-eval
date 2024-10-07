conda init
source ${7}
conda activate IGdetective

python3 ${6}/run_iterative_igdetective.py $1 $2
touch ${5}/igGene/${3}.${4}.txt
