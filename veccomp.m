function [all,dis]=veccomp(ij,n,F,X);   
for ji=1:n
               all(ji)=norm(F(ij,:)-F(ji,:));
               dis(ji)=norm(X(ij,:)-X(ji,:));
              end