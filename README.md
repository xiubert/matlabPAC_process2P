# matlabPAC_process2P

- ensure compiled code is added to MATLAB search path:  `compileRequiredCode.m` 

0.  Collect ad hoc FRA map during experiment:   `adhocFRAmap.m`
1.  Process raw data from individual animals:   `processAnimal2P.m`
2.  Compile processed data from multiple animals:   `compileCohortData.m`
3.  Plot compiled data from cohort: `plotCohortData.m`