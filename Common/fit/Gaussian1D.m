% A 1D, normalized, Gaussian function
%     
% FCN =  P(1) + P(3)/(P(2) * sqrt(2*pi)) * exp(- (x-P(4))^2 / (2*P(2)^2) )
%
% P(1) : background
% P(2) : sigma
% P(3) : area under the curve
% P(4) : mean
%
classdef Gaussian1D < FitFunction
    properties
        X; % array of x values
        guess; % the initial guess -> [bg sigma area mean]
        floatParams; % the float parameter array, eg. [1 0 1 1] -> float all parameters during the fit except for sigma
        nParams = 4; % the number of parameters in the function
    end
    methods
        
        function setX(self, X)
        % function setX(X)
        % 
        % Set the X values, an Nx1 array of x values
            if ~iscolumn(X)
                X = X';
            end
            setX@FitFunction(self, X);
        end
        
        function Y = eval(self, P)
        % function Y = eval(P)
        %
        % Evaluate the function using the parameter values in P
        %
        % Inputs
        % ------
        % P : 1D array
        %  parameter array -> [bg sigma area mean]
        %
        % Returns
        % -------
        % the function evaluated at each x value in the X array
        %
            assert(~isempty(self.X), 'Empty X array. Nothing to evaluate.\n');
            assert(length(P) == self.nParams, 'Parameter array does not have the proper length of %d, the length is %d.\n', self.nParams, length(P));
            Y = P(1)*ones(length(self.X), 1);
            Y = Y + P(1) + P(3)/(P(2) * sqrt(2*pi)) * exp(- (self.X - P(4)).^2 / (2*P(2)^2) );
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
            t2 = sqrt(2.0);
            t5 = 1.0 / sqrt(pi);
            t6 = P(2)^2;
            t12 = 1.0 / t6;
            t25 = t6^2;
            t10 = self.X - P(4);
            t11 = t10.^2;
            t15 = exp(-t11 * (t12 * 0.5));
            t16 = t2 * t5 / P(2) * t15 * 0.5;
            i = 1;
            j = 1;
            while (i <= self.nParams)
                if (i == 1) && self.floatParams(i)
                    jac(:,j) = 1.0;
                    j = j + 1;
                elseif (i == 2) && self.floatParams(i)                    
                    jac(:,j) = P(3) * t2 * t5 * 0.5 * t15 .* (-t12 + t11 ./ t25 );
                    j = j + 1;
                elseif (i == 3) && self.floatParams(i)
                    jac(:,j) = t16;
                    j = j + 1;
                elseif (i == 4) && self.floatParams(i)
                    jac(:,j) = P(3) * t16 / t6 .* t10;
                    j = j + 1;
                end
                i = i + 1;
            end
        end
        
        function newP = checkParams(~, bestParams)
        % function newP = checkParams(bestParams)
        %
        % Used during the LevMar fitting routine.
        %
        % Ensure that the angle, sigmax and aspect ratio values remain reasonable
        % for every fitting iteration.
        % 
        % Inputs
        % ------
        % bestParams : 1D array
        %  the current best-fit parameters
        %
        % Returns
        % -------
        % newP : 1D array
        %   the modifed paramter values
        %
            newP = bestParams;
        end
        
    end
    
end