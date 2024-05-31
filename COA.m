
function seeds=COA(A,seedsnumber,r,c)

degree=sum(A);
if c==0
    c=2*max(degree);
end
[Q,D]=eigs(A,r,'la');
seeds=[];
N=size(A,1);


for k=1:seedsnumber
    k
    bestconnectivity=1000000;
    bestindex=-1;
    for j=1:N
        if ismember(j,seeds)
            continue;
        end
        
        absu=sqrt(1-sum(Q(j,:).^2));
        R=[eye(r) Q(j,:)';zeros(1,r) absu];
        Lambda=[D zeros(r,1);zeros(1,r) -c];
        Z=R*Lambda*(R');
        [Q1,D1]=eig(full(Z));
        
        lambda=max(diag(D1));
        if lambda<bestconnectivity
            bestconnectivity=lambda;
            bestindex=j;
        end  
    end
    
    if bestindex==-1
        continue
    end
    j=bestindex;
    absu=sqrt(1-sum(Q(j,:).^2));
    R=[eye(r) Q(j,:)';zeros(1,r) absu];
    Lambda=[D zeros(r,1);zeros(1,r) -c];
    Z=R*Lambda*(R');
    [V1,D1]=eigs(Z,r+1,'lr');
    v=zeros(N,1);
    v(j)=1;
    u=v-sum(Q*diag(Q(j,:)),2);
    Q=[Q u]*V1;
    Q=Q(:,1:r);
    D=D1(1:r,1:r);    
    seeds=[seeds bestindex];
    
    
end

end








