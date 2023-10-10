# first time creating the venv

module load python gcc r/4.2.2 #best to do in home directory
virtualenv --no-download venv #you can call it whatever you want, here we call it venv
source venv/bin/activate

pip install --no-index --upgrade pip
pip install keras
pip install tensorflow

deactivate
