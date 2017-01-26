function [result]=ACMK(X,A,s,alpha,beta,gamma)
[m,n,nm]=size(A);
Z=eye(n);
c=length(unique(s));
e=1/12*ones(12,1);

for i=1:500
     K=zeros(n);
    for j=1:12
 K=K+e(j)*A(:,:,j);
    end
    
    
       Zold=Z;

     Z= (Z+Z')/2;
    
    D = diag(sum(Z));
    L = D-Z;
   
    [F, temp, ev]=eig1(L, c, 0);   
    
    parfor ij=1:n
              [all,dis]=veccomp(ij,n,F,X);
H=2*alpha*eye(n)+2*K;
H=(H+H')/2;
ff=beta*dis'+gamma/2*all'-2*K(:,ij);
              [Z(:,ij),err,lm] = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));   
          end
                   
          
     h=zeros(12,1);
for j=1:12
    h(j)=trace(A(:,:,j)-2*A(:,:,j)*Z+Z'*A(:,:,j)*Z);
end
 for j=1:12
 e(j)=(h(j)*sum(1./h))^(-2);    
 end
         
if i>10 &((norm(Z-Zold)/norm(Zold))<1e-3)
    
       break
end

end
actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
[result] = ClusteringMeasure( actual_ids,s);