%% Title : Linear least-squares regression
% Description :
% This file contains the code of the linear least-squares
% regression and the ridge regression. 

clear; close all; clc

%% Linear least-squares regression
% Generate the data: degree 1 
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

% Generate the data: degree 3 
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

% Linear least squares regression : example 1
tic
xc = [ones(n+1, 1), x'];
w = inv(xc' * xc) * xc' * ywn';
y_hat = xc*w;
toc 

for i = 1:n
    MSE1 = (ywn(i)-y_hat(i))^2;
end
MSE1 = 1/n * MSE1

figure
plot(x,y_hat, 'k', x,ywn, 'r*')
xlabel('x')
ylabel('y')
title('Regression line of data set 1')
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.90 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('y = %f * x + %f', w(2), w(1));
text(xt, yt, caption, 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');

% Linear least squares regression : example 2
tic
xc2 = [ones(n2+1, 1), x2', (x2.^2)', (x2.^3)'];
w2 = inv(xc2' * xc2) * xc2' * ywn2';
y_hat2 = xc2*w2;
toc

for i = 1:n
    MSE2 = (ywn2(i)-y_hat2(i))^2;
end
MSE2 = 1/n * MSE2

figure
plot(x2,y_hat2, 'k', x2,ywn2, 'r*')
xlabel('x')
ylabel('y')
title('Regression line of data set 2')
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.90 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('y = %f * x^3 +  %f * x^2 +  %f * x + %f', w2(4), w2(3), w2(2), w2(1));
text(xt, yt, caption, 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');

%% Ridge regression
% Ridge regression : with regularisation, example 1, lambda = 5 
tic
lambda = 5;
L = lambda.*eye(size(xc,2));
L(1) = 0;
w12 = inv(xc' * xc + L) * (xc') * ywn';
y_hat12 = xc*w12;
toc
figure
plot(x,y_hat12, 'k', x,ywn, 'r*')
xlabel('x')
ylabel('y')
title('Regression line of data set 1 with regularisation')
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.90 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('y = %f * x + %f', w12(2), w12(1));
text(xt, yt, caption, 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');

for i = 1:n
    MSE3 = (ywn(i)-y_hat12(i))^2;
end
MSE3 = 1/n * MSE3

% Ridge regression : with regularisation, example 2, lambda = 5 
tic
lambda = 5;
L2 = lambda.*eye(size(xc2,2));
L2(1) = 0;
w22 = inv(xc2' * xc2 + L2) * (xc2') * ywn2';
y_hat22 = xc2*w22;
toc
figure
plot(x2,y_hat22, 'k', x2,ywn2, 'r*')
xlabel('x')
ylabel('y')
title('Regression line of data set 2 with regularisation')
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.90 * (yl(2)-yl(1)) + yl(1);
caption = sprintf('y = %f * x^3 +  %f * x^2 + %f * x + %f', w22(4), w22(3) ,w22(2), w22(1));
text(xt, yt, caption, 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');

for i = 1:n
    MSE4 = (ywn2(i)-y_hat22(i))^2;
end
MSE4 = 1/n * MSE4