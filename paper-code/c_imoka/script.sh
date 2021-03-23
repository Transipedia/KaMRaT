#!/bin/bash
## To do for each CV

aggregated_matrix="./aggregated.matrix"  ## The final output of iMOKA computed from the training set
test_matrix="./test.json" ## the intial test matrix ( output of iMOKA_core create )

awk 'NR > 2 {print $1}' $aggregated_matrix > ./kmers.ls
singularity exec iMOKA_exended iMOKA_core extract -i ./kmers.ls -o ./test.matrix -s $test_matrix
### those two steps has to be done only once for each CV

## A model you want to test 
RF_model="./1_RF.pickle" 

singularity exec iMOKA_exended predict.py ./test.matrix $RF_model ./prediction_test.txt



