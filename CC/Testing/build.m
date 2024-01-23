function[L]=build(A)
     %hand coded cholesky of laplacian

     [n,m]=size(A);
     L=eye(n);
     for i=1:n-1
      L(i,i) = sqrt((i+1)/i);
      L(i,i+1) =  -1/L(i,i);
     end
     L(n,n)= sqrt(8);
     L(1,n)= -1/L(1,1);
     L(2,n) = -1/sqrt(6);
     for i=2:n-2
        L(i+1,n) = L(i,n)* sqrt((i)/(i+2));
     end
     L(n-1,n) = -sqrt(1+ 1/(n-1));

    L=sparse(L); 
