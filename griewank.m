function varargout = griewank(X)
% Griewank funcion 
%
%   GRIEWANK([x1, x2]) returns the value of the value of the Griewank
%   function at the specified points. [x1] and [x2] may be vectors.
%   The search domain is
%
%               -100 < x_i < 100
%
%   The global minimum is 
%
%               f(x1, x2) = f(0, 0) = 0.    
    
% Author: Rody P.S. Oldenhuis
% Delft University of Technology
% E-mail: oldenhuis@dds.nl
% Last edited 20/Jul/2009

    % if no input is given, return dimensions, bounds and minimum
    if (nargin == 0)
        varargout{1} = 2;  % # dims
        varargout{2} = [-100, -100]; % LB
        varargout{3} = [+100, +100]; % UB
        varargout{4} = [0, 0]; % solution
        varargout{5} = 0; % function value at solution
        
    % otherwise, output function value
    else
        
        % keep values in the search domain
        X(X < -100) = inf;  X(X > 100) = inf;
        
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
			sonuc2(i,1) = cos(X(i)/sqrt(i));
		end
		sonuc3=1;
		for i=1:D,
		  sonuc3=sonuc3*sonuc2(i,1);
		end
		  
        varargout{1} = sum(sonuc1)/4000 - sonuc3 + 1;
        
    end
     
end