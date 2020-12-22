% A 1D polynomial function of arbitrary order
%
% FCN = P(1) + P(2)*x + P(3)*y + P(4)*x^2 + P(5)*x*y + P(6)*y^2 + ... 
% 
% P(1) : background
% P(2) : x offset
% P(3) : 1st order coefficient
% P(4) : 2nd order coefficient
% ...
% P(n+2) : n'th order coefficient
%       
% WARNING! Be careful if you set P(2) ~= 0 and you have, say, P(3) and
% P(4) ~= 0. In this case, the P(4) factor will also have a linear term that 
% adds/subtracts with the contribution from the linear P(3) factor. Therefore,
% if P(2) ~= 0 and you want to use a polynomial of order n (coefficient P(n+2))
% then it is recommended that you should set the parameters P(3), P(4), ..., 
% P(n+1) to be equal to zero.
%
classdef Polynomial2D < FitFunction
    properties
        X; % array of x values
        guess; % the initial guess
        floatParams; % the float parameter array, eg. [1 0 1] -> float parameters 1, 3 during the fit; fix parameter 2
        nParams; % the number of parameters in the function
        order; % the order of the polynomial
    end
    properties(Hidden)
        npts; % the number of data points
    end
    methods
        
        function self = Polynomial2D(n)
        % Create a polynomial of a particular order
        %
        % Inputs
        % ------
        % n : integer
        %   the order of the polynomial
            self.order = n;
            if n == 0
                self.nParams = 1;
            else
                self.nParams = ((n+1)*(n+2))/2;
            end
            self.setFloatParams(ones(self.nParams, 1));
        end
        
        function setX(self, X)
        % function setX(X)
        % 
        % Set the X values.
        %
        % Inputs
        % ------
        % X : 2D matrix
        %  an Nx2 array of [x1 x2] coordinates (the independent variable)
            [r,c] = size(X);
            assert((r ~= 1) || (c ~= 1), 'X must be a Nx2 array, size=%dx%d', size(X))
            assert((r == 2) || (c == 2), 'X must be a Nx2 array, size=%dx%d', size(X))
            if c ~= 2
                X = transpose(X);
            end
            self.X = X;
            self.npts = size(X,1);
        end
        
        function y = eval(self, p)
        % function Y = eval(P)
        %
        % Evaluate the function using the parameter values in P.
        %
        % Inputs
        % ------
        % P : 1D array
        %   the parameter array
        %
        % Returns
        % -------
        % the function evaluated at each X value
        %
            assert(~isempty(self.X), 'Empty X array. Nothing to evaluate.\n');
            assert(length(p)==self.nParams, 'The length of the parameter array is incorrect. %d parameters are specified, require %d', length(p), self.nParams)            
            y = p(1)*ones(self.npts, 1);
            if self.nParams > 1
                for i=1:self.npts
                    index = 2;
                    for j=1:self.order
                        k = 0;
                        while k <= j
                            %fprintf('p(%d) * x^%d * y^%d\n', index, j-k, k)
                            y(i) = y(i) + p(index)*self.X(i,1)^(j-k)*self.X(i,2)^k;
                            index = index + 1;
                            k = k + 1;
                        end
                    end
                end
            end
        end
        
        function jac = jacobian(self, ~, jac)
        % function jac = jacobian(~, jac)
        %
        % The Jacobian matrix
        %
        % Inputs
        % ------
        % ~ : not used
        % 
        % jac : 2D array
        %  a pre-allocated matrix to contain the jacobian values
        %
        % Returns
        % -------
        % the Jacobian matrix
        %
            ifloat = 1;
            % P(1) contribution
            if self.floatParams(1)
                jac(:,ifloat) = 1.0;
                ifloat = ifloat + 1;
            end
            if self.nParams > 1
                idx = 2;
                for j=1:self.order
                    k = 0;
                    while k <= j
                        if self.floatParams(idx)
                            jac(:,ifloat) = self.X(:,1).^(j-k) .* self.X(:,2).^k;
                            ifloat = ifloat + 1;
                        end
                        idx = idx + 1;
                        k = k + 1;
                    end
                end
            end
        end        
        
        function s = toString(self)
        % function s = toString()
        %
        % Returns a string representation of this 2D polynomial
        % eg. if order=1 then returns 'P(1) + P(2)x + P(3)y'
            s = 'P(1)';
            index = 2;
            for j=1:self.order
                k = 0;
                while k <= j
                    s = strcat(s, sprintf(' + P(%d)', index));
                    if (j-k == 1)
                        s = strcat(s, 'x');
                    elseif (j-k > 1)
                        s = strcat(s, sprintf('x^%d', j-k));
                    end
                    if (k == 1)
                        s = strcat(s, 'y');
                    elseif (k > 1)
                        s = strcat(s, sprintf('y^%d', k));
                    end
                    index = index + 1;
                    k = k + 1;
                end
            end
        end
        
    end
end
