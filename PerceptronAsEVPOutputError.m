% Title : The perceptron with backpropagation
% Description :
% This code contains the solution of the perceptron with the
% Macaulay matrix method. 
% See also : 
% - Function_dataGeneration
% - Function_CreateSystem
% - constructmacaulayfast
% - SolveEVP
% - Function_plot_result
% - Function_plot_CostPerceptron1D
% - Function_Cost_Perceptron

clear; close all; clc
tic
% Initialise
[x,y] = Function_dataGeneration(1);
n = size(x,1);
nn = 2 * size(x,1) + size(x,2) + 1;

% Change the degree of the pad� function
% For a Degree 3
tx_num = [-14400 -1680 -24; 1 3 5];
tx_den = [14400 6480 264 1; 0 2 4 6];
% For a Degree 2
%tx_num = [-144 -12; 1 3];
%tx_den = [144 60 1; 0 2 4]; 
u = [ones(n,1), x];

P = Function_CreateSystem(u,y,tx_num,tx_den);

d = 7; % Degree 
K = SolveEVP(P,d,nn);
%Z = mxfnull(P,d);
%Mfast = constructmacaulayfast(P,d);

% Only keep the real solutions 
i = 0;
K_new = [];
for k = 1:size(K,2)
    K_r = K(:,k);
    column = K_r(~imag(K_r));
    if size(column,1) == size(K_r,1)
        i = i+1;
        K_new(:,i) = K_r;
    end 
end

count = 1;
for l = 1:size(K_new,2)
    for i = 1:n
        xx(l,i) = K_new(2,l)*u(i,2)+K_new(3,l);
    end
    
    pxk(l,:) = 48*xx(l,:).^5+3360*xx(l,:).^3+28800*xx(l,:);
    qxk(l,:) = 2*xx(l,:).^6+528*xx(l,:).^4+12960*xx(l,:).^2+28800;
    tx(l,:) = (pxk(l,:))./(qxk(l,:));
    
    MSE(l) = 0;
    MSE_equation(l) = 0;
    for ll = 1:size(pxk,2)
        MSE(l) = MSE(l) + (y(ll)-tx(l,ll))^2;
    end

    MSE(l) = MSE(l)/(size(u,1)); 
 
    % Print the cost 
%     fprintf('\nThe cost of result ')
%     disp(l)
%     fprintf('= ')
%     disp(MSE(l))
    
    count = count+1;
end
toc

% % Plot all the cost functions
% count = 1
% [X,Y,cost] = Function_plot_CostPerceptron1D_optimization(x,y);
% for l = 1:size(K_new,2)
%      % Plot the cost graph + results 
%     W1_new = [K_new(3,l), K_new(2,l)];
%     t = figure
%     mesh(X,Y,cost)
%     view([135 90 45])
%     xlabel('Weight')
%     ylabel('Bias')
%     zlabel('Cost')
%     hold on
%     plot3(W1_new(1,2),W1_new(1,1), MSE(l), 'r*')
%     hold off 
%     
%     saveas(t, sprintf('costfunction%d.png',count));
%     
%     count = count + 1;
% end

% Get the information  of best result 
imag_number = size(K,2)-size(K_new,2)
real_number = size(K_new,2)
best_MSE = min(MSE)
ind = find(MSE == min(MSE));
l = ind;
count = 1; 
Function_plot_result(1,x,y,[K_new(3,l), K_new(2,l)],count);

X = load('X.mat');
Y = load('Y.mat');
cost = load('Cost.mat');
W1_new = [K_new(3,l), K_new(2,l)]
t = figure
[X,Y,cost] = Function_plot_CostPerceptron1D_optimization(x,y);
mesh(X,Y,cost)
view([135 90 45])
xlabel('Weight')
ylabel('Bias')
zlabel('Cost')
hold on
plot3(W1_new(1,2),W1_new(1,1), MSE(l), 'r*')
hold off 

%saveas(t, sprintf('costfunctionda7%d.png',count));