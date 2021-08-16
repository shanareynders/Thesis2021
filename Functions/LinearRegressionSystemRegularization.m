function [P] = LinearRegressionSystemRegularization(x, y, lambda)
% [P] = LINEARREGRESSIONSYSTEM(x,y) returns the system that
    % can be used to create the Macaulay matrix for the ridge regression. 
    %
    %   x :  Input.
    %   y :  Output.
    %   lambda : regression parameter
    
numCells = size(x,2);
P = cell(numCells,1);
numColls = 1 + size(x,2);

% For all derivatives in w
for o = 1:size(x,2)
    temp = zeros(numColls, numColls);
    pos = 1;
    for i = 1:size(x,1)
       temp(pos,1) = (-2/size(x,1))*y(i)*x(i,o);
       pos = pos+1;
       for j = 1 : size(x,2)
            temp(pos,1) = x(i,j) * (2/size(x,1)) * x(i,o);
            temp(pos,1+j) = 1;
            pos = pos + 1; 
       end
    end
    if o ~= 1
        temp(pos,1) = 2*lambda/size(x,1);
        temp(pos,1+o) = 1;
    end
    P{o} = temp;
end 
end