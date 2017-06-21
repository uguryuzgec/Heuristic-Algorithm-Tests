function varargout = cube(X)
% Extended cube function
%
%   CUBE([x1, x2, ..., xn]) returns the value of the extended
%   cube function at the specified points. All [xi] may be 
%   vectors. 
%
%   The global minimum is 
%
%               f(x1, x2, ..., xn) = f(1, 1, ..., 1) = 0.

% Author: Rody P.S. Oldenhuis
% Delft University of Technology
% E-mail: oldenhuis@dds.nl
% Last edited 20/Jul/2009

    % if no input is given, return dimensions, bounds and minimum
    if (nargin == 0)
        varargout{1} = 2;  % # dims
        varargout{2} = [-10,-10]; % LB
        varargout{3} = [+10,+10]; % UB
        varargout{4} = [1,1]; % solution
        varargout{5} = 0; % function value at solution (too complicated)

    % otherwise, output function value
    else
        
        % keep values in the serach interval
        X(X < -100) = inf;     X(X > 100) = inf;
        
        % NOTE: orientation can not automatically be determined.
        % Defuaults to column sums....
        
        % split input vector X into X1, X2
% % %         X1 = X(1:end-1, :);    X2 = X(2:end, :);
       if size(X, 1) == 2
            x1 = X(1, :);        x2 = X(2, :);
        else
            x1 = X(:, 1);        x2 = X(:, 2);
        end
        
        % output columnsum
        varargout{1} = sum(100*(x2 - x1.^3).^2 + (1 - x1).^2, 1);
        
    end
     
end