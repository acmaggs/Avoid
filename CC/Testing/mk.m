function M=mk(N)
%periodic laplacian operator, with pinning on one site, to remove zero ev
M = diag(2*ones(1,N)) + diag(-1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1);
M(N,N)=10;
M(1,N)=-1;
M(N,1)=-1;
