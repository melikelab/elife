% VorDivergence.m
% script to get KL-Divergence scores with uniform model

% get uniform distribution of KL-Divergence
% likelihood = log ab^(c/a)/Gamma(c/a) + (c-1)1/n sum^n_{i=1} log(x_i) -
% b/n sum^n_{i=1} log (x_i^a)

% simulation parameters
Nsamples = 50;
Npoints = 200;
ReducedArea = 1;
Dim = ReducedArea*Npoints;

% statistical quantities
Areas = [];
% generate points and extract areas over and over again
for ii = 1:Nsamples
    xy = sqrt(Dim)*rand(Npoints,2);
    % make the boundaries periodic
    xy = [xy; [xy(:,1)+sqrt(Dim), xy(:,2)]; [xy(:,1)-sqrt(Dim), xy(:,2)]; ...
        [xy(:,1),xy(:,2)+sqrt(Dim)]; [xy(:,1),xy(:,2)-sqrt(Dim)]];
    DT = delaunayTriangulation(xy);
    [V,C] = voronoiDiagram(DT);
    VorDat = {V, C};
     
    % Calculate areas
    X = zeros(Npoints,1);
    for pt=1:Npoints
        xt = V(C{pt},1);
        yt = V(C{pt},2);
        X(pt) = abs(sum( (xt([2:end 1]) - xt).*(yt([2:end 1]) + yt))*0.5);
    end
    Areas = [Areas; X(~isnan(X))];
end

% fit theoretical distribution
LLfunc = @(a,b,c,X) a*b^(c/a) / gamma(c/a) * X.^(c-1) .* exp(-b*X.^a);

% determine coefficients
a = 1.07950;
b = 3.03226;
c = 3.31122;

% Entropy of the generalized Gaussian distribution
EntropyVal = log(b^(-1/a)*gamma(c/a)/a) + c/a + (1/a - c/a)*psi(c/a);

xx = 1e-2:0.01:10;
FuncVal = LLfunc(a,b,c,xx/ReducedArea);

% overlay histogram with expected fit
histogram(log(Areas),100,'Normalization','pdf','FaceColor','y','EdgeColor','w');
% hold on
% plot(xx,FuncVal/ReducedArea,'LineWidth',3);
% hold off

%% Voronoi plot
figure;
voronoi(xy(:,1),xy(:,2))
[Vv,Cc] = voronoin(xy, {'Qbb'});
color = {'r' 'b' 'g' 'm' 'c' 'y' 'k' 'w'} ;
for i = 1:length(Cc)
   fill(Vv(Cc{i},1),Vv(Cc{i},2),char(randsample(color,1))) ;
end

% define negative log likelihood function for losses
nLLfunc = @(a,b,c,X) -(c/a)*log(b) - log(a) + gammaln(c/a) - (c-1)*log(X) + b*X.^a;

Loss = nLLfunc(a,b,c,Areas/ReducedArea);

