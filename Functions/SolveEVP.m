function [K] = SolveEVP(P, d, n)
% [P] = SOLVEEVP(P,d,n) return the Vandermonde basis of the null space of
% the Macaulay matrix containing all solutions. 

    %   P:      Input system.
    %   d:      degree of the Macaulay matrix
    %   n:      Number of unknown variables/equations in the system

% STEP 1 : Construct Macaulay matrix
%Mfast = constructmacaulayfast(P,d)

% STEP 2 : Construct numerical basis of the null space of M via the SVD
%Z = null(Mfast);
Z = mxfnull(P,d);

% STEP 3 : Determine number of affine solutions and construct an affine basis
ma = 0;
m=1;

for k = 1:d
    nm = calculatenumberofmonomials(k,n)*m;
    W = Z(1:nm,:);
    r = rank(W,10e-10); % Only uses eigenvalues.
    
    if ma == r && ma > 0
        disp(['Number of affine solutions: ' num2str(ma)])
        break
    else
        ma = r;
    end
end
    
nm1 = calculatenumberofmonomials(k-1,n);
nm2 = calculatenumberofmonomials(k,n); 

% STEP 4 : Column compression
if size(Z,2) == ma 
    W = Z(1:nm2*m,:);
else
    [~,~,Q] = svd(Z(1:nm2*m,:));
    V = Z*Q;
    W = V(1:nm2*m,1:ma);
end

% STEP 5 : Construct the shift structure 
% Init the shift matrices:
%gx = [randn(n,1) eye(n)]; % randn vervangen door ones 
gx = [ones(n,1) eye(n)];
S1 = [eye(nm1*m) zeros(nm1*m,(nm2-nm1)*m)];
S2 = zeros(nm1*m,nm2*m);

% Construct the shift matrices 
K = createmonomialvectors(d-1,n);
for k = 1:nm1
    D = zeros(n,n+1); 
    for l = 1:n
        D(l,:) = [0 K(k,:)] + gx(l,:);
    end
    p = position(D(:,2:end));
    for l = 1:n
        S2((k-1)*m+1:k*m,(p(l)-1)*m+1:p(l)*m) = D(l,1)*eye(m);
    end
end

% Solve the EVP
[T,D] = eig(pinv(S1*W)*(S2*W));

% Retreive the Vandermonde basis K
K = W*T;
[rownumk,colnumk]=size(K);
for x = 1 : colnumk
    K(:,x) = K(:,x)/K(1,x);
end
end

