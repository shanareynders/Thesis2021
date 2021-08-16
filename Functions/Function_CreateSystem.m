function [P] = Function_CreateSystem(x,y, tx_num, tx_num)
% [P] = FUNCTION_CREATESYSTEM(x,y, tx_num, tx_den) returns the system that
    % can be used to create the Macaulay matrix for the output error. 
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

nCells = 2*N + nWeigths;
P = cell(nCells,1);
nCols = 2*N + nWeigths + 1;

% Degrees
nnum = max(tx_num(2,:));
nden = max(tx_den(2,:));
maxi = max(nnum, nden);

% All derivatives in tk 
d = 0;
for o = 1:N
    for t = 1:maxi
        s{t} = Function_Product(x(o,:), t, nCols);
    end
    
    temp_n = 2;
    for a = 1:size(tx_num,2)
        if tx_num(2,a) == 0
            temp_n = temp_n + 1;       
        else 
            temp_n = temp_n + size(s{tx_num(2,a)},1);
        end
    end
    temp = zeros(temp_n, nCols);
    
    % FILL IN
    temp(1,1) = 2*y(o)/N;
    
    temp(2,1) = -2/N;
    temp(2,N+(nWeigths+1)+o) = 1;
    pos = 3;
    
    for b = 1:size(tx_den,2)
        if tx_den(2,b) == 0
            temp(pos, 1) = tx_den(1,b);
            temp(pos, nWeigths+o+1) = 1;
            pos = pos+1;
        else
            ss = s{tx_den(2,b)};
            ss(:,1) = ss(:,1)*tx_den(1,b);
            ss(:,nWeigths+o+1) = 1;
            temp(pos:pos-1+size(ss,1),:) = ss;
            pos = pos+size(ss,1);
        end
    end
    
    d = d+1;
    P{d} = temp;
end

% All derivatives in lambda 
for o = 1:N
    for t = 1:maxi
        s{t} = Function_Product(x(o,:), t, nCols);
    end
    
    temp_n = 1;
    for a = 1:size(tx_num,2)
        if tx_num(2,a) ~= 0
            temp_n = temp_n + size(s{tx_num(2,a)},1);
        end
    end
    for a = 1:size(tx_den,2)
        if tx_den(2,a) ~= 0
            temp_n = temp_n + size(s{tx_den(2,a)},1);
        end
    end
    temp = zeros(temp_n, nCols);
    
    pos = 1;
    for b = 1:size(tx_den,2)
        if tx_den(2,b) == 0
            temp(pos, 1) = tx_den(1,b);
            temp(pos, nWeigths+N+o+1) = 1;
            pos = pos+1;
            
        else
            ss = s{tx_den(2,b)};
            ss(:,1) = ss(:,1)*tx_den(1,b);
            ss(:,nWeigths+N+o+1) = 1;
            temp(pos:pos-1+size(ss,1),:) = ss;
            pos = pos+size(ss,1);
        end
    end
    
    for b = 1:size(tx_num,2)
        if tx_num(2,b) == 0
            temp(pos, 1) = tx_num(1,b);
            pos = pos+1;
            
        else
            ss = s{tx_num(2,b)};
            ss(:,1) = ss(:,1)*tx_num(1,b);
            temp(pos:pos-1+size(ss,1),:) = ss;
            pos = pos+size(ss,1);
        end
    end
       
    d = d+1;
    P{d} = temp;
end

% All derivatives in the weigths
for t = 1:nWeigths 
    for o = 1:N
        for tt = 1:(maxi-1)
            s{tt} = Function_Product(x(o,:), tt, nCols);
        end
        
        temp_n = 1;
        for a = 1:size(tx_num,2)
            if tx_num(2,a)-1 > 0
                temp_n = temp_n + size(s{tx_num(2,a)-1},1);
            end
        end
        for a = 1:size(tx_den,2)
            if tx_den(2,a)-1 > 0
                temp_n = temp_n + size(s{tx_den(2,a)-1},1);
            end
        end
        temp = zeros(temp_n, nCols);
               
        pos = 1;
        for b = 1:size(tx_den,2)
            if tx_den(2,b)-1 == 0
                temp(pos, 1) = tx_den(1,b)*x(o,t);
                temp(pos, nWeigths+o+1) = 1;
                temp(pos, nWeigths+N+o+1) = 1;
                pos = pos+1;

            elseif tx_den(2,b)-1 > 0
                ss = s{tx_den(2,b)-1};
                ss(:,1) = ss(:,1)*tx_den(1,b)*tx_den(2,b)*x(o,t);
                ss(:,nWeigths+o+1) = 1;
                ss(:,nWeigths+N+o+1) = 1;
                temp(pos:pos-1+size(ss,1),:) = ss;
                pos = pos+size(ss,1);
            end
        end
        
        for b = 1:size(tx_num,2)
            if tx_num(2,b)-1 == 0
                temp(pos, 1) = tx_num(1,b)*x(o,t);
                temp(pos, nWeigths+o+1) = 1;
                pos = pos+1;

            elseif tx_num(2,b)-1 > 0
                ss = s{tx_num(2,b)-1};
                ss(:,1) = ss(:,1)*tx_num(1,b)*tx_num(2,b)*x(o,t);
                ss(:,nWeigths+o+1) = 1;
                temp(pos:pos-1+size(ss,1),:) = ss;
                pos = pos+size(ss,1);
            end
        end
                
        if o == 1
            nTemp = temp; 
        else
            nTemp = [nTemp;temp];
        end
    end
    d = d+1;
    P{d} = nTemp;
end
end

