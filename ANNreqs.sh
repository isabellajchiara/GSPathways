
salloc --time=1:00:00 --mem=30G --cpus-per-task=8 --account=def-haricots

module load gcc python r/4.2.2

source venv/bin/activate

cd projects

cd def-haricots

cd ich

cd MultipleCycles-main

R

