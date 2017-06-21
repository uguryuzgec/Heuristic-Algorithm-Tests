function varargout = sinenvsin(X)
% Sine envelope sine function
%
%   SINENVSIN([x1, x2, ..., xn]) returns the value of the sine envelope 
%   sine function at the specified points. [x1] and [x2] may be vectors. 
%
%   This is a generalized form of the Schaffer function--this function
%   results if the length of the input vector is 2. The search domain is
%
%               -100 < x_i < 100
%
%   The global minimum is 
%
%               f(x1, x2, ..., xn) = f(0, 0, ..., 0) = 0.

% Author: Rody P.S. Oldenhuis
% Delft University of Technology
% E-mail: oldenhuis@dds.nl
% Last edited 20/Jul/2009

    % if no input is given, return dimensions, bounds and minimum
    if (nargin == 0)
        varargout{1} = 2;  % # dims
        varargout{2} = [-20, -20]; % LB
        varargout{3} = [+20, +20]; % UB
        varargout{4} = [0, 0]; % solution
        varargout{5} = -1.4915; % function value at solution
    % otherwise, output function value
    else
        
        % keep all values in the search domain
        X(X < -20) = inf;  X(X > 20) = inf;
        
        % NOTE: orientation can not be determined automatically.
        % Defualts to column sums...
        % X1 = X(1:end-1, :);        X2 = X(2:end, :);
        
        % compute pythagorean sum only once
        % X12X22 = X1.^2 + X2.^2;
        
        % output columnsum
		% output function value
        D=size(X,2);
		for i=1:D-1,
			sonuc1(i,1) = (sin(sqrt(X(i+1)^2 + X(i)^2)- 0.5)^2)/((1 + 0.001*(X(i+1)^2 + X(i)^2))^2);
		end
        varargout{1} = (-1)*(sum(sonuc1)+0.5*(D-1));
    end
         
end