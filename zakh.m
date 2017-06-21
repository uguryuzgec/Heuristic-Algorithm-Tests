function varargout = zakh(X)
% 
% Zakharov function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.

  % if no input is given, return dimensions, bounds and minimum
    if (nargin == 0)
        varargout{1} = 2;        % # dims
        varargout{2} = [-5, -5]; % LB
        varargout{3} = [+10, +10]; % UB
        varargout{4} = [0, 0]; % solution
        varargout{5} = 0; % function value at solution
 % otherwise, output function value
    else
        % keep all values in the search domain
        X(X < -5) = inf;        X(X > 10) = inf;
        
        % split input vector X into x1, x2
        if size(X, 1) == 2
            x1 = X(1, :);        x2 = X(2, :);
        else
            x1 = X(:, 1);        x2 = X(:, 2);
        end
        
        % output function value
        s1 = 0;
        s2 = 0;
        D=size(X,2);
        for j = 1:D;
          s1 = s1+X(j)^2;
          s2 = s2+0.5*j*X(j);
        end
        varargout{1} = s1+s2^2+s2^4;;
    end
    
end