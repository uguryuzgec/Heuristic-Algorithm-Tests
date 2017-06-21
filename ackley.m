function varargout = ackley(X)
% Ackley funcion 
%
%   ACKLEY([x1, x2]) returns the value of the value of the Ackley
%   function at the specified points. [x1] and [x2] may be vectors.
%   The search domain is
%
%               -35 < x_i < 35
%
%   The global minimum is 
%
%               f(x1, x2) = f(0, 0) = 0

% Author: Rody P.S. Oldenhuis
% Delft University of Technology
% E-mail: oldenhuis@dds.nl
% Last edited 20/Jul/2009

    % if no input is given, return dimensions, bounds and minimum
    if (nargin == 0)
        varargout{1} = 2;  % # dims
        varargout{2} = [-35, -35]; % LB
        varargout{3} = [+35, +35]; % UB
        varargout{4} = [0, 0]; % solution
        varargout{5} = 0; % function value at solution

    % otherwise, output function value
    else
        
        % Keep all values in the search domain
        X(X < -35) = inf;   X(X > 35) = inf;
        
        % split input vector X into x1, x2
        if size(X, 1) == 2
            x1 = X(1, :);        x2 = X(2, :);
        else
            x1 = X(:, 1);        x2 = X(:, 2);
        end
        
        % output function value
        D=size(X,2);
		for i=1:D,
			sonuc1(i,1) = X(i)^2;
			sonuc2(i,1) = cos(2*pi*X(i));
		end
        varargout{1} = 20*(1 - exp(-0.2*sqrt(sum(sonuc1)/D)))...
				- exp((sum(sonuc2))/D) + exp(1);
    end
    
end