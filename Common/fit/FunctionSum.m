classdef FunctionSum < FitFunction
    properties        
        X; % the X values
        nParams; % the number of parameters in the function
        guess; % the initial guess
        floatParams;
        fcns; % array of functions to sum
    end
    properties(Hidden)
        npts; % the number of data points in X
        nFunctions; % the number of fit functions to sum
    end
    methods
        function self = FunctionSum(fcns)
        % function object = FunctionSum(fcns)
        %
        % Input
        % fcns : a cell array of FitFunctions
        %   eg. fcns = {Polynomial2D Gaussian2D}
        %
            assert(iscell(fcns), 'The functions must be in a cell array');
            self.nParams = 0;
            self.floatParams = [];
            self.guess = [];
            self.nFunctions = length(fcns);
            for i=1:self.nFunctions
                assert(isa(fcns{i}, 'FitFunction'), '%s is not a subclass of FitFunction', class(fcns{i}));
                self.nParams = self.nParams + fcns{i}.nParams;
                if i==1
                    size1 = size(fcns{i}.X);
                    if ~isempty(fcns{i}.X)
                        self.npts = size(fcns{i}.X, 1);
                        self.X = fcns{i}.X;
                    end
                else
                    assert(size1(1) == size(fcns{i}.X,1) && size1(2) == size(fcns{i}.X,2), 'The X array in %s does not have the same size as %s', class(fcns{i}), class(fcns{1}));
                end
                self.floatParams = [self.floatParams; fcns{i}.floatParams];                
                self.guess = [self.guess; fcns{i}.guess];
            end
            self.fcns = fcns;
            self.numFloating = sum(self.floatParams);         
        end
        
        function y = eval(self, P)
        % function y = eval(params)
        %
        % Evaluate the function using the parameter values in P
            assert(~isempty(self.X), 'Empty X array. Nothing to evaluate.\n');
            y = zeros(self.npts,1);
            idx = 1;
            for i=1:self.nFunctions
                y = y + self.fcns{i}.eval( P(idx: idx + self.fcns{i}.nParams - 1) );
                idx = idx + self.fcns{i}.nParams;
            end            
        end
        
        function setX(self, X)
        % function setX(X)
        %
        % Set the X values (the independent variable). The values must be
        % in column vectors.
            for i=1:self.nFunctions
                self.fcns{i}.setX(X);
            end
            self.npts = size(self.fcns{1}.X, 1);
            self.X = self.fcns{1}.X;
        end
        
        function setFloatParams(self, floatParams)
        % function setFloatParams(floatParams)
        %
        % Set the float-parameter array, eg. [1 1 1 0 1 1] -> float all
        % parameters during the fit except for parameter 4
            assert(length(floatParams)==self.nParams, 'The length of the float-parameter array must be %d, but the length is %d', self.nParams, length(floatParams));
            idx = 1;
            for i=1:self.nFunctions
                self.fcns{i}.setFloatParams(floatParams(idx: idx + self.fcns{i}.nParams - 1));
                idx = idx + self.fcns{i}.nParams;
            end
            setFloatParams@FitFunction(self, floatParams);
        end
        
        function setGuess(self, guess)
        % function setGuess(guess)
        %
        % Set the initial guess
            assert(length(guess)==self.nParams, 'The length of the guess array must be %d, but the length is %d', self.nParams, length(guess));
            idx = 1;
            for i=1:self.nFunctions
                self.fcns{i}.setGuess(guess(idx: idx + self.fcns{i}.nParams - 1));
                idx = idx + self.fcns{i}.nParams;
            end
            setGuess@FitFunction(self, guess);
        end
        
        function newP = checkParams(self, bestParams)
        % function out = checkParams(bestParams)
        % 
        % For example, to ensure that an angle is between -pi and pi
        %
        % Input
        % bestParams : the best parameter array
        %
            newP = bestParams;
            idx = 1;
            for i=1:self.nFunctions
                newP(idx: idx + self.fcns{i}.nParams - 1) = self.fcns{i}.checkParams(bestParams(idx: idx + self.fcns{i}.nParams - 1));
                idx = idx + self.fcns{i}.nParams;
            end
        end
        
        function jac = jacobian(self, P, jac)
        % function jac = jacobian(P, jac)
        % 
        % The Jacobian matrix
        %
        % Input
        % P : the parameter array
        % jac : pre-allocated 2D array
        %
            ij = 1;
            ip = 1;
            for i=1:self.nFunctions
                if self.fcns{i}.numFloating > 0
                    jac(:, ij: ij + self.fcns{i}.numFloating - 1) = self.fcns{i}.jacobian( P(ip: ip + self.fcns{i}.nParams - 1), jac(:, ij: ij + self.fcns{i}.numFloating - 1));
                end
                ij = ij + self.fcns{i}.numFloating;
                ip = ip + self.fcns{i}.nParams;
            end
        end
        
    end
end