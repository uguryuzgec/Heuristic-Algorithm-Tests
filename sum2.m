function varargout = sum2(X)
% 
% Sum Squares function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 20.
  % if no input is given, return dimensions, bounds and minimum
    if (nargin == 0)
        varargout{1} = 2;        % # dims
        varargout{2} = [-10, -10]; % LB
        varargout{3} = [+10, +10]; % UB
        varargout{4} = [0, 0]; % solution
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
        s = 0;
        D=size(X,2);
        for j = 1:D;
          s = s+j*X(j)^2; 
        end
        varargout{1} = s;
    end
    
end 