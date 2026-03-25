%%combine data from multiple animals 
load("C:\Users\JIC402\OneDrive - University of Pittsburgh\Data\AA0046\BPN1\AA0046_anmlROI_BPNTable.mat")
T1=GroupedTbl;
load("C:\Users\JIC402\OneDrive - University of Pittsburgh\Data\AA0046\BPN2\AA0046_anmlROI_BPNTable.mat")
T2=GroupedTbl;
load("C:\Users\JIC402\OneDrive - University of Pittsburgh\Data\AA0046\BPN3\AA0046_anmlROI_BPNTable.mat")
T3=GroupedTbl;

combinedTable = [T1; T2;T3];
GroupedTbl=combinedTable ;
