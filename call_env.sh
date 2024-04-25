

### subsequent sessions after your first run using create_env.sh
module load python gcc r/4.2.2
source venv/bin/activate

salloc --time=1:30:00 --mem=200G --cpus-per-task=4
cd projects
cd def-haricots
cd ich
cd MultipleCycles-Run
R
source("RunSims.R")



