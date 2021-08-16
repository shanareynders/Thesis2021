function [] = Function_plot_CostPerceptron1D(inputs, y, W1_new)
% [] = FUNCTION_PLOT_COSTPERCEPTRON1D(inputs, y, W1_new) makes a plot of
% the costfunction of a problem. 

    %   inputs:     the input data.
    %   y:          output data.
    %   W1_new:     The parameters of a certain solution that is plotted 
    %               together with the cost function.

    i = 1;
    m=size(y,1);
    x=[ones(m,1),inputs];
    
    for W = -10:0.01:10
        j = 1;
        for W2 = -10:0.01:10
            W1 = [W, W2];
            a1=x*W1' ;
            Z1 = Function_tanh(x*W1');

            %cost computations
            lambda = 2;
            [JJ, grad]= Function_Cost_Perceptron(W1,Z1,x,y);
            cost(i,j) = JJ;
            j = j+1;
        end
        i=i+1;
    end

    Z1 = Function_tanh(x*[W1_new(1,1), W1_new(1,2)]');
    [JJ, grad] = Function_Cost_Perceptron([W1_new(1,1), W1_new(1,2)],Z1,x, y);
    [X,Y] = meshgrid( -10:0.01:10, -10:0.01:10);
    
    figure()
    mesh(X,Y,cost)
    view([135 90 45])
    xlabel('Weight')
    ylabel('Bias')
    zlabel('Cost')
    hold on
    plot3(W1_new(1,2),W1_new(1,1), JJ, 'r*')
    hold off 
end

