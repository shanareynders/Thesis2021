% Title : The perceptron algorithm
% Description :
% This code contains the solution of the perceptron with the original learning algorithm 
% See also :
%   - Function_tanh
%   - Function_Cost_Perceptron
%   - Function_Cost_Perceptron
%   - Function_plot_result
%   - Function_plot_CostPerceptron1D

clear; close all; clc

[inputs, y] = Function_dataGeneration(1);

inputsize=size(inputs,2);
hiddensize=inputsize;
m=size(y,1);
epsilon_init = 0.16;

% Initializing weights % 
Length_in1  = inputsize;
Length_out1 = 1;

W1 = zeros(Length_out1,1+Length_in1);
W1 = rand(Length_out1, 1 + Length_in1)*2*epsilon_init - epsilon_init;
%W1 = [2,2]

% Forward computations % 
fprintf('\nFeedforward Using Neural Network ...\n')
x=[ones(m,1),inputs];
a1=x*W1' ;
Z1 = Function_tanh(x*W1');


% Cost computations
J = Function_Cost_Perceptron(W1,Function_tanh(x*W1'),x,y);

fprintf('\nThe cost is')
disp(J);

fprintf('\nTraining Neural Network... \n')
options = optimset('MaxIter',20,'Display','Iter-detailed');

costFunc = @(p) Function_Cost_Perceptron(p,Function_tanh(x*p'),x,y);

[W1_new, cost] = fminunc(costFunc, W1, options);

% Plot the results
figure()
Function_plot_result(1,inputs,y,W1_new,1);

% Plot the cost graph + results 
Function_plot_CostPerceptron1D(inputs,y,W1_new)