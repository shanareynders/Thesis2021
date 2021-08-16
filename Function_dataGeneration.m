function [x,y] = Function_dataGeneration(dimension)
% [x,y] = FUNCTION_DATAGENERATION(dimension) returns the input and output 
    % for a certain dimension together with a plot.
    %
    %   dimension : the dimension of the data

x = [3; -1; 2];
y = [1; -1; 1];

if dimension == 1
    labels=zeros(length(x),1);
    gscatter(x, labels,  y, 'br', '.',18);
    axis([min(x)-.1 max(x)+.1 -.05 .1])
    yticks(0)
    title('1-D Plot by the Class','FontSize',18)
elseif dimension == 2
    gscatter(x(:,1), x(:,2),  y, 'br', '.',18);
    axis([min(x(:,1))-.1 max(x(:,1))+.1 min(x(:,2))-.1 max(x(:,2))+.1])
    %yticks(0)
    title('2-D Plot by the Class','FontSize',18)
end
end

