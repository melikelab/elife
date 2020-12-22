% Abstract class for a FitFunction
%
classdef FitFunction < handle
    properties(Abstract)
        % the X values
        X;
        
        % the number of parameters in the function
        nParams; 

        % the initial guess
        guess;
        
        % 1D logical array indicating whether to float [0] or fix [1] 
        % each parameter during the fit. Example, floatParams = [0 1 0 1 1]
        % means to fix parameters 1, 3 and float parameteres 2, 4, 5
        floatParams;
    end
    properties(Hidden)
        numFloating; % the number of parameters that are floating
        eps = 1e-8; % the relative step size to use in the finite-difference estimates
    end
    methods
        
        function self = FitFunction()
        % function self = FitFunction()
        %
        % FitFunction construction. Initializes all floatParams = 1
            if ~isempty(self.nParams)
                self.setFloatParams(ones(self.nParams,1));
            end
        end

        function setX(self, X)
        % function setX(X)
        %
        % Set the X values (the independent variable). 
        %
        % Inputs
        % ------
        % X : 2D array of values
        %  The values must be in column vectors.        
            [r,c] = size(X);
            if c > r
                X = X(:);
            end        
            self.X = X;
        end

        function setGuess(self, guess)
        % function setGuess(guess)
        %
        % Set the initial guess.
        %
        % Inputs
        % ------
        % guess : 1D array
        %   the initial guess
        %
            assert(length(guess)==self.nParams, 'The length of the guess array must be %d, but the length is %d', self.nParams, length(guess));
            self.guess = guess(:);
        end

        function setFloatParams(self, floatParams)
        % function setFloatParams(floatParams)
        %
        % Set the float-parameter array.
        %
        % Inputs
        % ------
        % floatParams : 1D array of 1 and 0 (or true and false)
        %   eg. [1 1 1 0 1 1] -> float all parameters during the fit except for parameter 4
            assert(length(floatParams)==self.nParams, 'The length of the float-parameter array must be %d, but the length is %d', self.nParams, length(floatParams));
            self.floatParams = logical(floatParams(:));
            self.numFloating = sum(self.floatParams);
        end
        
        function out = checkParams(~, bestParams)
        % function out = checkParams(bestParams)
        % 
        % You must override this function if you want to modify the best fit
        % parameters for each iteration during the fit.
        % For example, to ensure that an angle is between -pi and pi
        %
        % Inputs
        % ------
        % bestParams : 1D array
        %   the best parameter array
            out = bestParams;
        end
        
        function setFiniteDifferenceFactor(self, eps)
        % function setFiniteDifferenceFactor(eps)
        %
        % Set the relative factor to increment each parameter in the
        % finite element difference estimation
        % 
        % Inputs
        % ------
        % eps : double
        %  the finite difference factor, epsilon
            self.eps = eps;
        end
        
        function jac = jacobian(self, params, jac)
        % function jac = jacobian(self, params, jac)
        %
        % You must override this function if you want to use analytical
        % formulas to calculate the Jacobian matrix. The default evaluation 
        % of the Jacobian is to use central finite element difference        
        %
        % Inputs
        % ------
        % params : 1D array
        %  the function parameters to use to evaluate the jacobian matrix
        %
        % jac : 2D array
        %  a pre-allocated array to contain the jacobian values
            n = length(params);
            dp = zeros(n,1);
            i = 1; ifloat = 1;
            while i <= n
                if self.floatParams(i)
                    h = self.eps * max(1, abs(params(i)));
                    dp(:) = 0;
                    dp(i) = 0.5*h;
                    jac(:,ifloat) = (self.eval(params + dp) - self.eval(params - dp))/h;
                    ifloat = ifloat + 1;
                end
                i = i + 1;
            end
        end        
    end
    
    methods(Abstract)
        % Evaluate the function using the parameter values in P
        y = eval(obj, P); 
    end
    
end