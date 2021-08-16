function [sigma] = Function_Product(x, degree, nCols)
% [sigma] = FUNCTION_PRODUCT(x, degree, nCols) calculates the product of a sum 
% with multiple dimensions. The degree can be chosen. 

    %   x:          Input data.
    %   degree:     The power of the sum
    %   nCols:      This determines the number of columns of the system in
    %               which this results needs to be implemented.

nWeigths = size(x,2);
n = nWeigths;

for i = 2:degree
    n = n*(nWeigths+i-1);
end
n = int8(n/factorial(degree)); 

sigma = zeros(n, nCols);

weigthVector(1, 1:degree) = ones(1,degree);
p = 2;
h = 1;

while weigthVector(p-1,degree) ~= nWeigths
    weigthVector(p, 1:degree) = weigthVector(p-1, 1:degree);
    if weigthVector(p,1) == nWeigths 
        weigthVector(p,2) = weigthVector(p,2) + 1;
        weigthVector(p,1) = weigthVector(p,2);
        for a = 2:degree
            if weigthVector(p,a) > nWeigths
                if a ~= degree
                    weigthVector(p,a+1) = weigthVector(p,a+1) + 1;
                    for b = 1:a
                        weigthVector(p,b) = weigthVector(p,a+1);
                    end
                end
            end
        end
    else
        weigthVector(p,1) = weigthVector(p,1)+1;
    end
    p = p +1;
end

o = 1;
for i = 1:size(weigthVector,1)
    constant = factorial(degree);
    
    h = histogram(weigthVector(i,:));
    values = h.Values;
    for k = 1:size(values, 2)
        constant = constant/factorial(values(k));
    end
     
    result = 1;
    for j = 1:degree
        result = result * x(weigthVector(i,j));
    end
    result = constant * result;
    
    sigma(o,1) = result;
    
    for s = 1:degree
        sigma(o,weigthVector(i,s)+1) = sigma(o,weigthVector(i,s)+1) + 1;
    end
    
    o = o+1;
end

end

