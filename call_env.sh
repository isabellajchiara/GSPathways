

### subsequent sessions after your first run using create_env.sh
salloc --time=1:30:00 --mem=90G --cpus-per-task=15
module load python gcc r/4.2.2
source venv/bin/activate
cd projects
cd def-haricots
cd ich
cd MultipleCycles-Run
R



