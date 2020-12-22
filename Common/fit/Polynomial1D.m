% A 1D polynomial function of arbitrary order
%
% FCN = P(1) + P(3)*(x-P(2)) + P(4)*(x-P(2))^2 + ... + P(n+2)*(x-P(2))^n
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
classdef Polynomial1D < FitFunction
    properties
        X; % array of x values
        guess; % the initial guess
        floatParams; % the float parameter array, eg. [1 0 1] -> float parameters 1, 3 during the fit; fix parameter 2
        nParams; % the number of parameters in the function
        order; % the order of the polynomial
    end
    properties(Hidden)
        xOffset; % allocate an array of (x-P(2)) values
    end
    methods
        
        function self = Polynomial1D(n)
        % function self = Polynomial1D(n)
        %
        % Create a polynomial of a particular order
        %
        % Inputs
        % 
        % n : integer
        %   the order of the polynomial
        %
            self.order = n;
            if n == 0
                self.nParams = 1;
                self.setFloatParams(1);
            else
                self.nParams = n + 2;
                fp = ones(self.nParams, 1);
                fp(2) = 0; % by default, don't float the x offset in the fit   
                self.setFloatParams(fp);
            end
        end

        function setX(self, X)
        % function setX(X)
        % 
        % Set the X values.
        %
        % Inputs
        % ------
        % X : 1D array
        %  the x values (the independent variable)
            self.xOffset = zeros(length(X),1);
            setX@FitFunction(self, X);
        end
        
        function Y = eval(self, P)
        % function Y = eval(P)
        %
        % Evaluate the function using the parameter values in P.
        %
        % Inputs
        % ------
        % P : 1D array
        %   the parameter array -> [background x0 c1 c2 c3 ...]
        %
        % Returns
        % -------
        % the function evaluated at each X value
        %
            assert(~isempty(self.X), 'Empty X array. Nothing to evaluate.\n');
            assert(length(P) == self.nParams, 'Parameter array does not have the proper length of %d, the length is %d.\n', self.nParams, length(P));
            Y = P(1)*ones(length(self.X),1);
            if self.nParams > 1
                self.xOffset = self.X - P(2);
                i = 2;
                while i < self.nParams
                    Y = Y + P(i+1)*self.xOffset.^(i-1);
                    i = i + 1;
                end                
            end
        end
        
        function setFloatParams(self, floatParams)
        % function setFloatParams(floatParams)
        %
        % Inputs
        % ------
        % floatParams : 1D array of 0's and 1's (false and true)
        %   Set which parameters are floating during the fit, [1 1 0 1] means
        %   float parameters 1, 2, 4 and fix parameter 3
        %
            %fp = logical(floatParams);
            %if (length(fp) > 1) && (fp(2) ~= 0) && any(fp(3:end-1))
            %    recommend = zeros(size(floatParams));
            %    recommend([1 2 end]) = 1;
            %    fprintf('WARNING! Poly.setFloatParams() :: There are "dangerous" parameters floating in the fit\nYou specified %s. Recommended is %s if you want P(2) floating\n', mat2str(floatParams), mat2str(recommend));
            %end
            setFloatParams@FitFunction(self, floatParams);
        end
        
        function jac = jacobian(self, P, jac)
        % function jac = jacobian(P, jac)
        %
        % The Jacobian matrix
        %
        % Inputs
        % ------
        % P : 1D array
        %  the function parameters to use to evaluate the jacobian matrix
        %
        % jac : 2D array
        %  a pre-allocated matrix to contain the jacobian values
        %
        % Returns
        % -------
        % the Jacobian matrix
        %
            ifloat = 1;
            % the P(1) contribution
            if self.floatParams(1)
                jac(:,ifloat) = 1.0;
                ifloat = ifloat + 1;
            end
            if self.nParams > 1
                self.xOffset = self.X - P(2);
                if self.floatParams(2)
                    ixo = self.floatParams(1) + 1;
                    ifloat = ifloat + 1;
                    jac(:,ixo) = 0;
                end
                i = 2; 
                while i < self.nParams
                    % the P(3), P(4), ... P(n+2) contribution
                    if self.floatParams(i+1)
                        jac(:,ifloat) = self.xOffset.^(i-1);
                        ifloat = ifloat + 1;
                    end
                    % the P(2) contributions
                    if self.floatParams(2) && (P(i+1) ~= 0)
                        jac(:,ixo) = jac(:,ixo) + (1.0-i)*P(i+1)*self.xOffset.^(i-2);
                    end
                    i = i + 1;
                end
            end
        end
        
        function s = toString(self)
        % function s = toString()
        %
        % Returns a string representation of this 1D polynomial
        % eg. if order=2 then returns 'P(1) + P(3)(x-P(2)) + P(4)(x-P(2))^2'
            s = 'P(1)';
            index = 3;
            for j=1:self.order
                s = strcat(s, sprintf(' + P(%d)(x-P(2))', index));
                if (j > 1)
                    s = strcat(s, sprintf('^%d', j));
                end
                index = index + 1;
            end            
        end
        
    end
    
end