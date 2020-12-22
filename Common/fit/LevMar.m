% Determine the function parameters that minimize the chi square by solving a 
% non-linear least-squares problem using the Levenberg–Marquardt algorithm.
%
% Solve the system of equations using Cholesky decomposition. If the 
% J'WJ + lambda*diag(J'WJ) matrix is not positive definite then the Cholesky
% decomposition will fail. If this is the case, then solve eq. 3.2 in 
%   "Levenberg--Marquardt algorithm: implementation and theory"
%   Jorge J. More, Conference on numerical analysis, Dundee, UK, 28 Jun 1977
%   http://www.osti.gov/scitech/biblio/7256021
% using QR decomposition
%
classdef LevMar < handle
    properties
        fcn; % the fit function
        bestParams; % the best fit parameters
        paramUncerts; % the uncertainty of the best fit parameters
        maxIter; % the maximum number of fitting iterations
        iter; % the number of fitting iterations
        reducedChisqr; % the reduced chisqr value
        tol; % tolerance
        covMatrix; % the covariance matrix
        weights;
        yValues;
        residuals; % the residuals, y - y(params)
    end
    properties(Hidden)
        lambda; % the scaling factor
        numFloating; % the number of parameters floating during the fit
        numFree; % the number of degrees of freedom
        numParams; % the number of parameters
        initGradNorm; % the initial value of the gradient Norm
        npts; % the number of data points
        jac; % the Jacobian matrix
        jwj; % evaluation of the J'WJ matrix multiplication, J=Jacobian, W=Weights, '=transpose
        gradient; % evaluation of the J'*W*residuals matrix multiplication, J=Jacobian, W=weights, '=transpose
        choleskyL; % the lower triangual matrix of the Cholskey Decomposition
        deltaP; % the amount to increment each parameter
        paramTol; % the parameter differences between fit iterations
        acceptStep;
        trialChisqr;
        rejectCount;
        trialParams;
        trialResiduals;
    end
    methods
        function self = LevMar(fcn, yValues, yUncerts)
        % function object = LevMar(fcn, yValues, yUncerts)
        %
        % Inputs
        % fcn : subclass of a FitFunction object
        %   A fcn brings the X data, the initial guess and whether a
        %   parameter if floating or fixed during the fit
        %
        % yValues : 1D array
        %   The Y values that the fit should converge to
        %
        % yUncerts (optional argument) : the uncertainty of each Y value
        %    Used to calculate the weights for the fit, weight=1/yUncert^2
        %    If not specified then the yUncert(:) = 1 so weights(:) = 1

            if ~iscolumn(yValues)
                self.yValues = yValues';
            else
                self.yValues = yValues;
            end

            assert(nargin > 1, 'You must specify a fit function and the Y values');
            assert(isa(fcn, 'FitFunction'), 'LevMar only accepts a subclass of a FitFuntion as the "fcn" input variable');
            assert(~isempty(self.yValues), 'The yValues array is empty');
            assert(size(self.yValues,2) == 1, 'The yValues must be a Nx1 array. Got a %dx%d array', size(self.yValues))
            assert(size(fcn.X,1) == length(self.yValues), 'The X and yValues arrays are not the same length, %d ~= %d', size(fcn.X,1), length(self.yValues))
            
            self.npts = length(yValues);
            if nargin == 2                
                self.weights = ones(self.npts, 1);
            else
                assert(self.npts == length(yUncerts), 'The yValues and yUncerts arrays have different lengths, %d ~= %d', self.npts, length(yUncerts));
                self.weights = (1.0 ./ yUncerts.^2);
            end
            
            self.fcn = fcn;
        end
        
        function exit_code = solve(self, varargin)
        % Determine the parameters that minimize the chi square.
        %
        % Input
        % varargin : a key-value pair
        %   Allowed keys are 'maxIter', 'tol'
        %   If varargin is not specified then the default values are used
        %   example: solve('maxIter', 100, 'tol', 1e-6)
        %   If you misspell a key (case sensitive) then the specified value will be ignored
        %
        % Return
        % exit_code : int
        %   Value will be >0 if the fit converged or <0 if something went wrong
        %   exit_code =  3 , fit converged by the reduced chi-square condition
        %   exit_code =  2 , fit converged since the relative change in each parameter is small
        %   exit_code =  1 , fit converged by the norm(gradient) condition
        %   exit_code = -1 , maximum number of fitting iterations reached
        %   exit_code = -2 , NaN or Inf occurred
        %   exit_code = -3 , if the parameter increment is rejected a certain number of times in a row then stop fitting
        
            % initialize all matrices and values
            self.initialize(varargin);
            
            % get the initial chi square value
            self.residuals = self.getResiduals(self.bestParams);
            self.reducedChisqr = self.getChiSquare(self.residuals);
            
            exit_code = 0;
            A = [];
            self.rejectCount = 0;
            self.acceptStep = true;
            self.iter = 0;
            while ~exit_code
                
                % update the Levenberg-Marquardt equation
                if self.acceptStep                
                    self.updateLM();
                end
                
                % update the value of lambda
                self.updateLambda();
                
                % Solve the system of equations by Cholesky decomposition.
                % If the J'WJ + lambda*diag(J'WJ) matrix is not positive definite then the Cholesky
                % decomposition will fail. If this is the case then solve eq. 3.2 in 
                % "Levenberg--Marquardt algorithm: implementation and theory"
                % Jorge J. More, Conference on numerical analysis, Dundee, UK, 28 Jun 1977
                % http://www.osti.gov/scitech/biblio/7256021
                if ~self.solveCholesky()
                    if self.acceptStep || isempty(A)
                        A = [self.jac; sqrt(self.lambda)*eye(self.numFloating)];
                        B = [self.residuals; zeros(self.numFloating, 1)];
                    else
                        A(self.npts+1:end,:) = sqrt(self.lambda)*eye(self.numFloating);
                    end
                    
                    warningstate1 = warning('off','MATLAB:nearlySingularMatrix');
                    warningstate2 = warning('off','MATLAB:singularMatrix');
                    warningstate3 = warning('off','MATLAB:rankDeficientMatrix');
                    
                    % solve by QR decomposition
                    self.deltaP = A \ B;

                    warning(warningstate1)
                    warning(warningstate2)
                    warning(warningstate3)                    
                end
                
                % increment the parameters
                self.stepParameters();
                
                % check if the new parameters should be accepted or rejected
                self.acceptStep = self.checkStep();
                
                % if it is a good step then accepts 
                if self.acceptStep
                    self.bestParams = self.fcn.checkParams(self.trialParams);
                    self.residuals = self.trialResiduals;
                    self.reducedChisqr = self.trialChisqr;
                end
                
                % exit_code will be >0 if the fit converged or <0 if something went wrong
                exit_code = self.checkStatus();
                
                self.iter = self.iter + 1;
            end
            
            % calculate the covariance matrix and the best fit parameters uncertainties
            if (exit_code >= -1)
                if self.makeCovarianceMatrix();
                    self.getParamUncerts();
                else
                    self.paramUncerts = ones(size(self.bestParams)) * NaN;
                end
            end
        end
        
    end
    
    methods(Hidden)
        
        function initialize(self, varargin)
            
            % parse the varargin variables
            key   = {'maxIter',        'tol'};
            value = [      200,        1e-8];
            dict = containers.Map(key,value);
            for i=1:2:length(varargin{1})-1
                if isKey(dict, varargin{1}{i})
                    dict(varargin{1}{i}) = varargin{1}{i+1};
                end
            end            
            self.maxIter = dict('maxIter');
            self.tol = dict('tol');
            
            % numbers
            self.numFloating = self.fcn.numFloating;
            self.numParams = length(self.fcn.guess);
            assert(self.numParams > 0, 'The guess has not been set yet');
            self.numFree = self.npts - self.numFloating;
            assert(self.numFree > 0, 'The number of free parameters is <= 0 (i.e. there are more floating fitting variables than data points)');

            % matrices
            self.bestParams = self.fcn.guess;
            self.jac = zeros(self.npts, self.numFloating);
            self.jwj = zeros(self.numFloating);
            self.choleskyL = zeros(self.numFloating);
            self.paramTol = ones(self.numFloating, 1);
        end
        
        function updateLambda(self)
           if (self.iter == 0)
                % initialize lambda to be 0.1% of the largest J'WJ diagonal value
                maxDiagJWJ = -realmax;
                i = 1;
                while i <= self.numFloating
                    maxDiagJWJ = max(maxDiagJWJ, abs(self.jwj(i,i)));
                    i = i + 1;
                end
                self.lambda = max(0.001, maxDiagJWJ*0.001);
                % the initial gradient norm
                self.initGradNorm = max(abs(self.gradient));
           elseif self.acceptStep
               % the new parameters are good, decrease lambda
               self.lambda = max(self.lambda*0.1, 1e-10);
           else
               % the new parameters are good, increase lambda
               self.lambda = min(self.lambda*10.0, 1e10);
           end
        end
        
        function accept = checkStep(self)
        % specify the conditions that are necessary to accept a step
        
            accept = false;
            if self.trialChisqr < self.reducedChisqr
                accept = true;
                self.rejectCount = 0;

                % how much each parameter has changed
                i = 1; j = 1;
                while i <= self.numParams
                    if self.fcn.floatParams(i)
                        if self.bestParams(i) ~= 0                            
                            self.paramTol(j) = abs(self.deltaP(j)/self.bestParams(i));                            
                        else
                            self.paramTol(j) = realmax;
                        end
                        j = j + 1;
                    end
                    i = i + 1;
                end
                
            else
                self.rejectCount = self.rejectCount + 1;
            end
        end
        
        function stepParameters(self)
            
            % increment the parameters that are allowed to vary during the fit
            self.trialParams = self.bestParams;
            i = 1;
            j = 1;
            while i <= self.numParams 
                if self.fcn.floatParams(i)
                    self.trialParams(i) = self.trialParams(i) + self.deltaP(j);
                    j = j + 1;
                end
                i = i + 1;
            end
            
            % get the trial chisqr
            self.trialResiduals = self.getResiduals(self.trialParams);
            self.trialChisqr = self.getChiSquare(self.trialResiduals);

        end
        
        function updateLM(self)
        % Update LM equation: (J'WJ + lambda*diag(J'WJ)) = J'*W*residuals
            self.jac = self.fcn.jacobian(self.bestParams, self.jac);
            self.gradient = self.jac' * (self.weights .* self.residuals);
            i = 1;
            while i <= self.numFloating
                j = 1;
                while j <= i
                    self.jwj(i,j) = self.weights' * (self.jac(:,i) .* self.jac(:,j));
                    self.jwj(j,i) = self.jwj(i,j);
                    j = j + 1;
                end
                i = i + 1;
            end
        end
        
        function chisqr = getChiSquare(self, residuals)
            chisqr = (self.weights' * residuals.^2)/self.numFree;
        end
        
        function res = getResiduals(self, params)
            res = self.yValues - self.fcn.eval(params);
        end
        
        function success = makeCholeskyL(self, A)
        % make the lower-triangular matrix for Cholesky decomposition, A=LL'            
            self.choleskyL(:) = 0;
            i = 1;
            while i <= self.numFloating
                j = i;
                while j <= self.numFloating
                    summ = A(i,j);
                    k = 1;
                    while k < j
                        summ = summ - self.choleskyL(i,k)*self.choleskyL(j,k);
                        k = k + 1;
                    end
                    if (i == j)
                        if (summ <= 0.0) % then the matrix is not positive definite
                            success = false;
                            return
                        end
                        self.choleskyL(i,i) = sqrt(summ);
                    else
                        self.choleskyL(j,i) = summ/self.choleskyL(i,i);
                    end
                    j = j + 1;
                end
                i = i + 1;
            end
            success = true;
        end
        
        function x = solveLowerTriangluar(self, b)
        % Solve the system of equations, LL'x=b, where L=lower triangular
            
            % Solve L*y=b, storing y in x
            x = zeros(self.numFloating,1);
            i = 1;
            while i <= self.numFloating
                s = b(i);
                k = i;
                while k > 0
                    s = s - self.choleskyL(i,k) * x(k);
                    k = k - 1;
                end
                x(i) = s / self.choleskyL(i,i);
                i = i + 1;
            end
            
            % Solve L'*x=y
            i = self.numFloating;
            while i > 0
                s = x(i);
                k = i + 1;
                while k <= self.numFloating
                    s = s - self.choleskyL(k,i) * x(k);
                    k = k + 1;
                end
                x(i) = s / self.choleskyL(i,i);
                i = i - 1;
            end
        end
        
        function success = solveCholesky(self)
        % Solve the system of equations, Ax=b, using Cholesky Decomposition
            modJWJ = self.jwj;
            i = 1;
            while i <= self.numFloating
                modJWJ(i,i) = self.jwj(i,i)*(1.0 + self.lambda);
                i = i + 1;
            end
            
            success = self.makeCholeskyL(modJWJ);
            if success
                self.deltaP = self.solveLowerTriangluar(self.gradient);
            end
        end
        
        function success = makeCovarianceMatrix(self)
        % calculate the inverse of J'WJ            
            self.covMatrix = NaN * ones(self.numFloating);
            success = self.makeCholeskyL(self.jwj);            
            if success
                e = zeros(self.numFloating, 1);
                col = 1;
                while col <= self.numFloating
                    e(col) = 1;
                    self.covMatrix(:,col) = self.solveLowerTriangluar(e);
                    e(col) = 0;
                    col = col + 1;
                end
            end
        end
        
        function getParamUncerts(self)
        % determine the parameter uncertainties from the covariance matrix
            self.paramUncerts = zeros(self.numParams,1);
            i = 1;
            j = 1;
            while i <= self.numParams
                if self.fcn.floatParams(i)
                    self.paramUncerts(i) = sqrt(self.covMatrix(j,j)*self.reducedChisqr);
                    j = j + 1;
                end
                i = i + 1;
            end
        end
        
        function exit_code = checkStatus(self)
        % Test if the fit converged or if a problem occurred
        % exit_code < 0 , problem occurred, stop fitting
        % exit_code = 0 , continue fitting
        % exit_code > 0 , fit converged, stop fitting
            
            exit_code = 0;
            
            % check for NaN or Inf in the step size
            if any(isnan(self.deltaP)) || any(isinf(self.deltaP))
                exit_code = -2;
                return
            end
            
            % if the parameter increment is rejected too many times in a row then stop fitting
            if self.rejectCount > 20
                exit_code = -3;
                return
            end

            % check if the current norm is a factor of tolerance smaller than the initial norm
            if max(abs(self.gradient)) < self.tol * self.initGradNorm
                exit_code = 1;
                return
            end            
            
            % check if the relative change in each parameter is less than the tolerance
            if (self.iter > 2) && all(self.paramTol < self.tol)
                exit_code = 2;
                return
            end

            % check if the reduced chi square is less than tol
            if self.reducedChisqr < self.tol
                exit_code = 3; 
                return
            end                
            
            % check if the maximum number of iterations was exceeded
            if self.iter == self.maxIter -1
                exit_code = -1;
                return 
            end

        end
        
    end
end