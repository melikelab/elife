% A 2D, normalized, elliptical Gaussian function
%     
% FCN = P(1) + P(4)/(2*pi*P(7)*P(6)^2) * 
%       exp( -1/2*( (x-P(2))*cos(P(5)) - (y-P(3))*sin(P(5)) )^2/P(6)^2
%            -1/2*( (x-P(2))*sin(P(5)) + (y-P(3))*cos(P(5)) )^2/(P(7)*P(6))^2 )
%
% P(1) : background
% P(2) : central x value, x0
% P(3) : central y value, y0
% P(4) : area under the curve
% P(5) : rotation angle (in radians)
% P(6) : sigma value in the x direction
% P(7) : aspect ratio, sigmaY/sigmaX (or equivalently fwhmY/fwhmX)
%
classdef Gaussian2D < FitFunction
    properties
        X; % Nx2 array of [x,y] pixel coordinates
        guess; % the initial guess -> [bg x0 y0 area angle sigmax aspect]
        floatParams; % the float parameter array, eg. [1 1 1 1 0 1 1] -> float all parameters during the fit except for the angle
        nParams = 7; % the number of parameters in the function
    end
    methods
        
        function setX(self, X)
        % function setX(X)
        % 
        % Set the X values, an Nx2 array of [x,y] coordinates
            assert(size(X,2)==2, 'The X values must be in a Nx2 array. Received a %dx%d array', size(X));
            setX@FitFunction(self, X);
        end
        
        function X = generateX(~, x0, y0, nx, ny)
        % function generateX(x0, y0, nx, ny)
        %
        % Generate an array of Nx2 pixel coordinates. 
        % Useful if you just want to model data for a particular parameter array. 
        %
        % For example, if x0=85, y0=125, nx=1, ny=2 then X would represent
        % a 15x2 array with the following values,
        %  (84,123) (85,123) (86,123)
        %  (84,124) (85,124) (86,124)
        %  (84,125) (85,125) (86,125)
        %  (84,126) (85,126) (86,126)
        %  (84,127) (85,127) (86,127)
        %
        % Inputs
        % ------
        % x0 : the central pixel coordinate in the x dimension
        % y0 : the central pixel coordinate in the y dimension
        % nx : the number of pixels to use on both sides of x0, in the x dimension
        % ny : the number of pixels to use on both sides of y0, in the y dimension
        %
        % Returns
        % -------
        % the X array
        %
            idx = 1;
            X = zeros((2*nx+1)*(2*ny+1), 2);
            for y=-ny:ny
                for x=-nx:nx                
                    X(idx,1) = x0 + x;
                    X(idx,2) = y0 + y;
                    idx = idx + 1;
                end
            end
        end        
        
        function Y = eval(self, P)
        % function Y = eval(P)
        %
        % Evaluate the function using the parameter values in P
        %
        % Inputs
        % ------
        % P : 1D array
        %  parameter array -> [bg x0 y0 area angle sigmax aspect]
        %
        % Returns
        % -------
        % the function evaluated at each [x y] coordinate in the X array
        %
            assert(~isempty(self.X), 'Empty X array. Nothing to evaluate.\n');
            assert(length(P) == self.nParams, 'Parameter array does not have the proper length of %d, the length is %d.\n', self.nParams, length(P));
            dx = self.X(:,1) - P(2);
            dy = self.X(:,2) - P(3);
            Y = P(1)*ones(size(self.X(:,1), 1), 1);
            Y = Y + P(4) / (2.0 * pi * P(6) * (P(7)*P(6))) * ...
                exp(-0.5 * ( (dx * cos(P(5)) - dy * sin(P(5))).^2 / P(6)^2 ...
                           + (dx * sin(P(5)) + dy * cos(P(5))).^2 / (P(7)*P(6))^2 ) );
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
            cosine = cos(P(5));
            sine = sin(P(5));
            tc0 = 1.0 / P(6);
            tc1 = 1.0 / P(7);
            x0 = self.X(:,1) - P(2);
            y0 = self.X(:,2) - P(3);
            ta1 = x0 * cosine - y0 * sine;
            ta2 = x0 * sine + y0 * cosine;
            ta3 = ta2 .* ta2 * (tc1 * tc1);
            ta4 = (ta1 .* ta1 + ta3) * (tc0 * tc0);
            ta5 = (P(4) * tc0 * tc0 * tc1 / (2.0 * pi)) * exp(-0.5 * ta4 ) ;

            i = 1;
            j = 1;
            while (i <= self.nParams)
                if (i == 1) && self.floatParams(i)
                    jac(:,j) = 1.0;
                    j = j + 1;
                elseif (i == 2) && self.floatParams(i)
                    jac(:,j) = (ta5 * (tc0 * tc0)) .* (ta1 * cosine + ta2 * (tc1 * tc1 * sine));
                    j = j + 1;
                elseif (i == 3) && self.floatParams(i)
                    jac(:,j) = (ta5 * (tc0 * tc0)) .* (ta2 * (tc1 * tc1 * cosine) - ta1 * sine);
                    j = j + 1;
                elseif (i == 4) && self.floatParams(i)
                    jac(:,j) = ta5 * (1.0 / P(4));
                    j = j + 1;
                elseif (i == 5) && self.floatParams(i)
                    jac(:,j) = ((1.0 - tc1 * tc1) * tc0 * tc0) * (ta5 .* ta1 .* ta2);
                    j = j + 1;
                elseif (i == 6) && self.floatParams(i)
                    jac(:,j) = (ta5 * tc0) .* (ta4 - 2.0);
                    j = j + 1;
                elseif (i == 7) && self.floatParams(i)
                    jac(:,j) = (ta5 * tc1) .* ((tc0 * tc0) * ta3 - 1.0);
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
            % convert the angle to be between -pi and pi
            newP(5) = bestParams(5) - 2.0*pi*(round(bestParams(5)/(2.0*pi))); 
            % The value of sigmax and the aspect ratio depend on the value of the angle.
            % For example, the function evaluated at theta=pi/2, sigmax=2, aspect=0.5
            % is equal to the function evaluated at theta=0 (or pi, or -pi), sigmax=1, aspect=2
            % Force the angle to be between -pi/4 < theta < pi/4 and change the value of
            % sigmax and the aspect ration when the value of the angle is within certain ranges.
            if newP(5) >= 0.75*pi
                newP(5) = newP(5) - pi;
            elseif (newP(5) > 0.25*pi) && (newP(5) < 0.75*pi)
                newP(5) = newP(5) - 0.5*pi;
                newP(6) = newP(6) * newP(7);
                newP(7) = 1.0 / newP(7);
            elseif (newP(5) > -0.75*pi) && (newP(5) < -0.25*pi)
                newP(5) = newP(5) + 0.5*pi;
                newP(6) = newP(6) * newP(7);
                newP(7) = 1.0 / newP(7);
            elseif (newP(5) <= -0.75*pi)
                newP(5) = newP(5) + pi;
            end
        end
        
    end
    
end