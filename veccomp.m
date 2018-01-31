function [all,dis]=veccomp(ij,n,F,X);   
for ji=1:n
               all(ji)=(norm(F(ij,:)-F(ji,:)))^2;
               dis(ji)=(norm(X(ij,:)-X(ji,:)))^2;
              end
