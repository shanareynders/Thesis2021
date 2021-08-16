function [P] = Function_CreateSystemEquationError(x,y, tx_num, tx_den)
% [P] = FUNCTION_CREATESYSTEM(x,y, tx_num, tx_den) returns the system that
    % can be used to create the Macaulay matrix for the equation error. 
    %
    %   tx_num:  [coefficients; degrees] of the factors of the numerator 
    %   of the approximated tangent hyperbolic activation function
    %   tx_num:  [coefficients; degrees] of the factors of the denominator 
    %   of the approximated tangent hyperbolic activation function
    %   x :  Input.
    %   y :  Output.

% x = [b w1 w2 ... wP]
% Initialisations
nWeigths = size(x,2);
N = size(x,1);

nCells = nWeigths;
P = cell(nCells,1);
nCols = nWeigths + 1;

% Degrees
nnum = max(tx_num(2,:));
nden = max(tx_den(2,:));
%maxi = max(nnum, nden);
maxi = nnum + nden;

d = 0;
% All derivatives in the weigths
for t = 1:nWeigths 
    for o = 1:N
        for tt = 1:maxi
            s{tt} = Product(x(o,:), tt, nCols);
        end
        
        temp_n1 = size(tx_num,2)+size(tx_den,2);
        temp_n = 0;
             
        for a = 1:size(tx_num,2)
            if tx_num(2,a)-1 > 0
                for b = 1:size(tx_num,2)
                    temp_n = temp_n + size(s{tx_num(2,a)-1+tx_num(2,b)},1);
                end
                for b = 1:size(tx_den,2)
                    temp_n = temp_n + size(s{tx_num(2,a)-1+tx_den(2,b)},1);
                end
            elseif  tx_num(2,a)-1 == 0
                for b = 1:size(tx_num,2)
                    if tx_num(2,b) > 0
                        temp_n = temp_n + size(s{tx_num(2,b)},1);
                    else
                        temp_n = temp_n + 1;
                    end
                end
                for b = 1:size(tx_den,2)
                    if tx_den(2,b) > 0
                        temp_n = temp_n + size(s{tx_den(2,b)},1);
                    else
                        temp_n = temp_n + 1;
                    end
                end
            end
        end
        for a = 1:size(tx_den,2)
            if tx_den(2,a)-1 > 0
                for b = 1:size(tx_num,2)
                    temp_n = temp_n + size(s{tx_den(2,a)-1+tx_num(2,b)},1);
                end
                for b = 1:size(tx_den,2)
                    temp_n = temp_n + size(s{tx_den(2,a)-1+tx_den(2,b)},1);
                end
            elseif  tx_den(2,a)-1 == 0
                for b = 1:size(tx_num,2)
                    if tx_num(2,b) > 0
                        temp_n = temp_n + size(s{tx_num(2,b)},1);
                    else
                        temp_n = temp_n + 1;
                    end
                end
                for b = 1:size(tx_den,2)
                    if tx_den(2,b) > 0
                        temp_n = temp_n + size(s{tx_den(2,b)},1);
                    else
                        temp_n = temp_n + 1;
                    end
                end
            end
        end
        temp = zeros(temp_n, nCols);
               
        pos = 1;
        for b = 1:size(tx_den,2)
            if tx_den(2,b) == 0
                for c = 1:size(tx_den,2)
                    if tx_den(2,c)-1 == 0
                        temp(pos,1) = tx_den(1,b)*2*y(o)^2*x(o,t)*tx_den(1,c)*tx_den(2,c);
                        pos = pos+1;
                    elseif tx_den(2,c)-1 > 0
                        ss = s{tx_den(2,c)-1};
                        ss_new = ss;
                        ss_new(:,1) = ss_new(:,1)*tx_den(1,b)*2*y(o)^2*x(o,t)*tx_den(1,c)*tx_den(2,c);
                        temp(pos:pos-1+size(ss_new,1),:) = ss_new;
                        pos = pos+size(ss,1);
                    end
                end
                for c = 1:size(tx_num,2)
                    if tx_num(2,c)-1 == 0
                        temp(pos,1) = tx_den(1,b)*2*y(o)*x(o,t)*tx_num(1,c)*tx_num(2,c);
                        pos = pos+1;
                    elseif tx_num(2,c)-1 > 0
                        ss = s{tx_num(2,c)-1};
                        ss_new = ss;
                        ss_new(:,1) = ss_new(:,1)*tx_den(1,b)*2*y(o)*x(o,t)*tx_num(1,c)*tx_num(2,c);
                        temp(pos:pos-1+size(ss_new,1),:) = ss_new;
                        pos = pos+size(ss,1);
                    end
                end
            elseif tx_den(2,b) > 0
                for c = 1:size(tx_den,2)
                    if tx_den(2,c)-1 >= 0
                        ss = s{tx_den(2,b)+tx_den(2,c)-1};
                        ss_new = ss;
                        ss_new(:,1) = ss_new(:,1)*tx_den(1,b)*2*y(o)^2*x(o,t)*tx_den(1,c)*tx_den(2,c);
                        temp(pos:pos-1+size(ss_new,1),:) = ss_new;
                        pos = pos+size(ss,1);
                    end
                end
                for c = 1:size(tx_num,2)
                    if tx_num(2,c)-1 >= 0
                        ss = s{tx_den(2,b)+tx_num(2,c)-1};
                        ss_new = ss;
                        ss_new(:,1) = ss_new(:,1)*tx_den(1,b)*2*y(o)*x(o,t)*tx_num(1,c)*tx_num(2,c);
                        temp(pos:pos-1+size(ss_new,1),:) = ss_new;
                        pos = pos+size(ss,1);
                    end
                end
            end
        end
        
        for b = 1:size(tx_num,2)
            if tx_num(2,b) == 0
                for c = 1:size(tx_den,2)
                    if tx_den(2,c)-1 == 0
                        temp(pos,1) = tx_num(1,b)*2*y(o)*x(o,t)*tx_den(1,c)*tx_den(2,c);
                        pos = pos+1;
                    elseif tx_den(2,c)-1 > 0
                        ss = s{tx_den(2,c)-1};
                        ss_new = ss;
                        ss_new(:,1) = ss_new(:,1)*tx_num(1,b)*2*y(o)*x(o,t)*tx_den(1,c)*tx_den(2,c);
                        temp(pos:pos-1+size(ss_new,1),:) = ss_new;
                        pos = pos+size(ss,1);
                    end
                end
                for c = 1:size(tx_num,2)
                    if tx_num(2,c)-1 >= 0
                        temp(pos,1) = tx_num(1,b)*2*x(o,t)*tx_num(1,b)*tx_num(2,b);
                        pos = pos+1;
                    elseif tx_num(2,c)-1 > 0
                        ss = s{tx_num(2,c)-1};
                        ss_new = ss;
                        ss_new(:,1) = ss_new(:,1)*tx_num(1,b)*2*x(o,t)*tx_num(1,c)*tx_num(2,c);
                        temp(pos:pos-1+size(ss_new,1),:) = ss_new;
                        pos = pos+size(ss,1);   
                    end
                end
            elseif tx_num(2,b) > 0
                for c = 1:size(tx_den,2)
                    if tx_den(2,c)-1 >= 0
                        ss = s{tx_num(2,b)+tx_den(2,c)-1};
                        ss_new = ss;
                        ss_new(:,1) = ss_new(:,1)*tx_num(1,b)*2*x(o,t)*tx_den(1,b)*tx_den(2,b);
                        temp(pos:pos-1+size(ss_new,1),:) = ss_new;
                        pos = pos+size(ss,1);
                    end
                end
                for c = 1:size(tx_num,2)
                    if tx_num(2,c)-1 >= 0
                        ss = s{tx_num(2,b)+tx_num(2,c)-1};
                        ss_new = ss;
                        ss_new(:,1) = ss_new(:,1)*tx_num(1,b)*2*x(o,t)*tx_num(1,b)*tx_num(2,b);
                        temp(pos:pos-1+size(ss_new,1),:) = ss_new;
                        pos = pos+size(ss,1);
                    end
                end
            end
        end        
        
        if o == 1
            nTemp = temp; 
        else
            nTemp = [nTemp;temp];
        end
    end
    d = d+1;
    nTemp(:,1) = nTemp(:,1)/N;
    P{d} = nTemp;
end
end

