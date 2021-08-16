%% Title : Linear least-squares regression + Ridge regression with 
% Macaulay marix method 
% Description :
% This Matlab code solves the linear regression problem with and without
% regularization with the Macaulay matrix method
% See also functions : 
%   - LinearRegressionSystem
%   - LinearRegressionSystemRegularization
%   - SolveEVP

clear; close all; clc

%% Without regularization : Example 1
n = 100;
x = 0:0.1:10;
y = 0.5 * x + 0.9;

figure; plot(x,y)

noise = wgn(n+1,1,5);
ywn = y + noise.';

plot(x, ywn, '*r')
xlabel('x')
ylabel('y')
title('Data set 1')

xc = [ones(n+1, 1), x'];

tic
P = LinearRegressionSystem(xc, ywn);

d = 1;
K = SolveEVP(P,d,2);
toc
[rownumk,colnumk]=size(K);
for i = 1 : colnumk
    K(:,i) = K(:,i)/K(1,i);
end

figure; plot(x,x*K(3,1)+K(2,1), 'r', x,ywn, 'b*')

% Without regularization : Example 2
n2 = 100;
x2 = 0:0.1:10;
y2 = 0.8 .* x2.^3 + 0.6 .* x2.^2 + 0.5 .* x2 + 0.9;

figure; plot(x2,y2)

noise2 = wgn(n2+1,1,40);
ywn2 = y2 + noise2.';

plot(x2, ywn2, '*r')
xlabel('x')
ylabel('y')
title('Data set 2')

xc2 = [ones(n2+1, 1), x2', (x2.^2)', (x2.^3)'];
tic
P = LinearRegressionSystem(xc2, ywn2);

d = 1;
K = SolveEVP(P,d,4);
toc
[rownumk,colnumk]=size(K);
for i = 1 : colnumk
    K(:,i) = K(:,i)/K(1,i);
end

figure; plot(x2,x2.^3*K(5,1)+ x2.^2*K(4,1)+x2*K(3,1)+K(2,1), 'r', x2,ywn2, 'b*')

%% With regularization : Example 1
n = 100;
x = 0:0.1:10;
y = 0.5 * x + 0.9;

figure; plot(x,y)

noise = wgn(n+1,1,5);
ywn = y + noise.';

plot(x, ywn, '*r')
xlabel('x')
ylabel('y')
title('Data set 1')

xc = [ones(n+1, 1), x'];

tic
P = LinearRegressionSystemRegularization(xc,ywn,5);

d = 1;
K = SolveEVP(P,d,2);
toc
[rownumk,colnumk]=size(K);
for i = 1 : colnumk
    K(:,i) = K(:,i)/K(1,i);
end

figure; plot(x,x*K(3,1)+K(2,1), 'r', x,ywn, 'b*')

% Without regularization : Example 2
n2 = 100;
x2 = 0:0.1:10;
y2 = 0.8 .* x2.^3 + 0.6 .* x2.^2 + 0.5 .* x2 + 0.9;

figure; plot(x2,y2)

noise2 = wgn(n2+1,1,40);
ywn2 = y2 + noise2.';

plot(x2, ywn2, '*r')
xlabel('x')
ylabel('y')
title('Data set 2')

xc2 = [ones(n2+1, 1), x2', (x2.^2)', (x2.^3)'];
tic
P = LinearRegressionSystemRegularization(xc2,ywn2,5);

d = 1;
K = SolveEVP(P,d,4);
toc
[rownumk,colnumk]=size(K);
for i = 1 : colnumk
    K(:,i) = K(:,i)/K(1,i);
end

figure; plot(x2,x2.^3*K(5,1)+ x2.^2*K(4,1)+x2*K(3,1)+K(2,1), 'r', x2,ywn2, 'b*')
