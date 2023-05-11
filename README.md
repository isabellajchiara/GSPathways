# These scripts simulate three breeding cycles. They are used to investigate when to select parents for the next cycle in a GS scenario, and which generation to include in model training. 

To run replications of the simultion, use the "RunReplicationsMC.R" Script. This script will collect genetic values, model performance, and variance data.

Start with "CycleOne". This script calls CycleTwo and CycleTwo calls CycleThree.

Try different models by sourcing: random forest(RF), neural net(ANN), support vector machine(SVM) or ridge regression(RRBLUP)

Models can use random(RD) or stratified clusters(SC) training set. 

Parent selections can be made at F2 or F5 by sourcing "ParentSelectionsF2" or "ParentSelectionsF5".
Relocate the parent selection script to select parents at a different generation

You may replace genetic maps and haplotype data with your own genetic map and SNP data.

