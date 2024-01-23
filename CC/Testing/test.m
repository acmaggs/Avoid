clear
n=7
M=mk(n);
X=build(M);
dim=500000;
r = normrnd(0,1,[n,dim]);

s=X'*r;
%should be cose to original matrix
cor=s*s'/dim


r = normrnd(0,1,[n,dim]);
inv2=X \r ;
%inverse of M: correlation matrix
cor2= inv2 *inv2'/dim;

%should be close to (eye)
prod=(M*cor2 +cor2*M)/2




%test my solver, and compare to simple matlab solution random vector "r"
r= normrnd(0,1,[n,1]);
r=1:n;
r=sqrt(r');
d=squeeze(diag(X));
o=squeeze(diag(X,1));
c=X(:,n);
out(1:n)=0;

out(n) = r(n)/d(n);
out(n-1)= (r(n-1)-out(n)*o(n-1) )/ d(n-1);

for i = 0:n-3
out(n-2-i)= (r(n-2-i)-out(n-1-i)*o(n-2-i) - c(n-2-i)*out(n))/ d(n-2-i);
end
%my solution
my_solution=out
inv2= X  \r;
matlab_solution=inv2'
diagonal=d'
off_diagonal=o'
column=c'
