function h = Function_tanh(z)
% [h] = FUNCTION_TANH(z) Computes the tanh function

    %   z:  Weigthed sum of the input

h = (exp(z)-exp(-z))./(exp(z)+exp(-z)); 
end