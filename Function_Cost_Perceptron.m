function [J, grad] = Function_Cost_Perceptron(W1,Z1, x, y)
    % [J, grad] = FUNCTION_COST_PERCEPTRON(W1,Z1, x, y) calculates the cost
    % and the updated parameters 
    %
    %   W1:  The parameters
    %   Z1:  The output of the perceptron, the activation function of the
    %   weigthed sum
    %   x :  Input.
    %   y :  Output.

    % Calculation of the cost function
    mse = sum((Z1 - y).^2);
    J = mse/size(x,1);

    %gradient for back propogation
    c = 0.1;
    grad = W1 + (1/2)*c*(Z1-y).*x;
end

