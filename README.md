# These scripts run three breeding cycles. They are used to investigate when to select parents for the next cycle, and which generation to include in model training. 

Here you can run 3 breeding cycles in a row

Start with "CycleOne". This script calls CycleTwo and CycleTwo calls CycleThree.

Try different models by sourcing: random forest(RF), neural net(ANN), support vector machine(SVM) or ridge regression(RRBLUP)

Models can use random(RD) or stratified clusters(SC) training set. 

Parent selections can be made at F2 or F5 by sourcing "ParentSelectionsF2" or "ParentSelectionsF5".
Relocate the parent selection script to select parents at a different generation

You may replace genetic maps and haplotype data with your own genetic map and SNP data.
