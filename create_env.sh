
module load python gcc r/4.2.2
virtualenv --no-download venv
source venv/bin/activate

pip install --no-index --upgrade pip
pip install keras
pip install tensorflow

deactivate

### subsequent sessions
module load python gcc r/4.2.2
source ANN/bin/activate


