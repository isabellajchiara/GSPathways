These scripts simulate any number of breeding cycles under a genomic selection framework using the pedigree breeding method.

You can choose between three prediction models, RRBLUP, RF or SVM. For each model, you may choose which generation to train at, and whether or not you would like a randomly sampled training population or a training population sampled via stratified clustering. 

Here's how to use this tool:

1) Load your SNP data into your WD, or generate genotype data using AlphaSimR. For us the genotype file is called "haplotypesSNPS.RData".
2) Load your genetic map into your WD. For us, the map file is called "genMapSNPs.RData" and contains the location of each SNP in Morgans.
(view these files if you need an example for formatting your own data)
4) Source the file RunSims.R and follow the prompts
5) Usethe DataVisualization files to read the simulation outputs
   - Genetic values at each generation
   - Model performance (Accuracy) at each generation when the model is implemented
   - Genetic variance at each generation
   - BV and EBV at each generation when the model is implemented
   - Genotype (alleles coded 1 or 0) at each generation
     

