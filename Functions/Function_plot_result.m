function [] = Function_plot_result(dimension, x, y, W, k)
% [] = FUNCTION_PLOT_result(dimension, x, y, W, k) makes a plot 
% of the result of a problem. 

    %   dimension:  The input dimension
    %   x:          Input data.
    %   y:          Output data.
    %   W:          The weigths 
    %   k:          A number to save the image automatically

if dimension == 1
    labels=zeros(length(x),1);

    t = figure
    gscatter(x, labels,  y, 'br', '.',18);
    axis([min(x)-.1 max(x)+.1 -.05 .1])
    yticks(0)
    title('1-D Plot by the Class','FontSize',18)
    hold on
    plot(-W(1)/W(2), 0, 'g*');
    legend('y = -1','y = 1','boundary')
    hold off
    
    fig = sprintf('figda7%d',k);
    %saveas(t, sprintf('fig%d',k), 'png');
    saveas(t, fig, 'png');
end
end

