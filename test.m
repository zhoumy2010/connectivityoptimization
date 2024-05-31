load('matrix.mat')
seeds=COA(netscience,20,20,30);
[result]=evallargesteigenvalue(netscience,seeds,1);
plot(result(:,1),result(:,2),'-*')
