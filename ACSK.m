function [result]=ACSK(X,K,s,alpha,gamma)
[m,n]=size(K);
Z=eye(n);
c=length(unique(s));

for i=1:100
  Zold=Z;
     Z= (Z+Z')/2;
    D = diag(sum(Z));
    L = D-Z;
   
    [F, temp, ev]=eig1(L, c, 0);   
    
    
         parfor ij=1:n

              [all,dis]=veccomp(ij,n,F,X);

H=2*alpha*eye(n)+2*K;
H=(H+H')/2;
ff=dis'+gamma/2*all'-2*K(:,ij);
              [Z(:,ij),err,lm] = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));   

          end
if i>10 &((norm(Z-Zold)/norm(Zold))<1e-3)
       break
end

end

actual_ids= kmeans(F, c, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');
           
[result] = ClusteringMeasure( actual_ids,s);
 