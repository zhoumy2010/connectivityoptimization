function [result]=evallargesteigenvalue(A,seeds,step)
L=length(seeds);
result=[];
for i=1:step:L
    A(seeds(1:i),:)=0;
    A(:,seeds(1:i))=0;
    [V,D]=eigs(A);
    val=D(1,1);
    result=[result;i val];
    
end



end