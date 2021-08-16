% Title : The perceptron with backpropagation
% Description :
% This code contains the solution of the perceptron with the
% backpropagation learning algorithm. 
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
nn = size(x,2)+1;

% Change the degree of the padé function
tx_num = [-14400 -1680 -24; 1 3 5];
tx_den = [14400 6480 264 1; 0 2 4 6];
u = [ones(n,1), x];

P = Function_CreateSystemEquationError(u,y,tx_num,tx_den);

d = 25; % Degree 
%Mfast = constructmacaulayfast(P,d);
K = SolveEVP(P,d,nn);
%Z = mxfnull(P,d);

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
%     Function_plot_result(1,x,y,[K_new(3,l), K_new(2,l)],count);
    
    for i = 1:n
        xx(l,i) = K_new(3,l)*u(i,2)+K_new(2,l);
        
    end
    
    pxk(l,:) = 48*xx(l,:).^5+3360*xx(l,:).^3+28800*xx(l,:);
    qxk(l,:) = 2*xx(l,:).^6+528*xx(l,:).^4+12960*xx(l,:).^2+28800;
    tx(l,:) = (pxk(l,:))./(qxk(l,:));
    
    MSE(l) = 0;
    MSE_equation(l) = 0;
    for ll = 1:size(tx,2)
        MSE(l) = MSE(l) + (y(ll)-tx(l,ll))^2;
    end

    MSE(l) = MSE(l)/(size(u,2)); 
    
    % Print the cost 
%     fprintf('\nThe cost of result ')
%     disp(l)
%     fprintf('= ')
%     disp(MSE(l))
    
    count = count+1;
end
toc

% % Plot the cost function
% count = 1
% [X,Y,cost] = Function_plot_CostPerceptron1D_optimization(x,y);
% for l = 1:size(K_new,2)
%      % Plot the cost graph + results 
%      W1_new = [K_new(3,l), K_new(2,l)];
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

% Get the information of the best result
imag_number = size(K,2)-size(K_new,2)
real_number = size(K_new,2)
best_MSE = min(MSE)
ind = find(MSE == min(MSE));
l = ind;
count = 1; 
Function_plot_result(1,x,y,[K_new(2,l), K_new(3,l)],count);

W1_new = [K_new(2,l), K_new(3,l)]
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
% 
% saveas(t, sprintf('costfunctiond21%d.png',count));