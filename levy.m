function varargout = levy(X)
% 
% Levy function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n =2.
 if (nargin == 0)
        varargout{1} = 2;        % # dims
        varargout{2} = [-10, -10]; % LB
        varargout{3} = [+10, +10]; % UB
        varargout{4} = [1, 1]; % solution
        varargout{5} = 0; % function value at solution
 % otherwise, output function value
    else
        % keep all values in the search domain
        X(X < -10) = inf;        X(X > 10) = inf;
        
        % split input vector X into x1, x2
        if size(X, 1) == 2
            x1 = X(1, :);        x2 = X(2, :);
        else
            x1 = X(:, 1);        x2 = X(:, 2);
        end
        
        % output function value
        D=size(X,2);
        for i = 1:D,
            z(i) = 1+(X(i)-1)/4;
        end
        s = sin(pi*z(1))^2;
        for i = 1:D-1,
            s = s+(z(i)-1)^2*(1+10*(sin(pi*z(i)+1))^2);
        end
        varargout{1} = s+(z(D)-1)^2*(1+(sin(2*pi*z(D)))^2);
    end
    
end 