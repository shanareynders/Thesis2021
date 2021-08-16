function [X,Y,cost] = Function_plot_CostPerceptron1D_optimization(inputs, y)
% [] = FUNCTION_PLOT_COSTPERCEPTRON1D_optimization(inputs, y) makes a plot 
% of the costfunction of a problem. 

    %   inputs:     the input data.
    %   y:          output data.

    i = 1;
    m=size(y,1);
    x=[ones(m,1),inputs];
    
    for W = -10:0.01:10
        j = 1;
        for W2 = -10:0.01:10
            W1 = [W, W2];
            a1=x*W1' ;

            pxk = 48*a1.^5+3360*a1.^3+28800*a1;
            qxk = 2*a1.^6+528*a1.^4+12960*a1.^2+28800;
            tx = (pxk)./(qxk);

            %cost computations
            lambda = 2;
            [JJ, grad]= Function_Cost_Perceptron(W1,tx,x,y);
            cost(i,j) = JJ;
            j = j+1;
        end
        i=i+1;
    end
    
    [X,Y] = meshgrid( -10:0.01:10, -10:0.01:10);
    t = figure
    mesh(X,Y,cost)
    view([180 0])
    xlabel('Weight')
    ylabel('Bias')
    zlabel('Cost')
    hold on
    plot3(W1_new(1,2),W1_new(1,1), cost_solution, 'r*')
    hold off 
%     
%     saveas(t, sprintf('costfunction%d.png',k));
end

