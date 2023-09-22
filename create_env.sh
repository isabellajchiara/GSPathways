virtualenv --no-download tensorflow
source tensorflow/bin/activate

pip install --no-index --upgrade pip
pip install keras
pip install tensorflow

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

R



