% HybridKLGammaTest.m
% script to test the hybrid KL divergence estimator against a generalized
% gamma distribution

% define parameters
a = 1.07950;
b = 3.03226;
c = 3.31122;
UniCoeff=[a,b,c];
%N = 20050;

% log PDF of generalized gamma distribution in log area space
PDFgamlog = @(p,y) p(3)/p(1) * log(p(2)) + log(p(1)) - gammaln(p(3)/p(1)) ...
    + p(3)*y - p(2)*exp(y*p(1));
UNIgamlog = @(y) PDFgamlog(uniCoeff,y);

% use voronoi simulation to get the gamma stretched histogram
% simulation parameters
Nsamples = 100;
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

% Entropy of the generalized Gaussian distribution
EntropyVal = log(b^(-1/a)*gamma(c/a)/a) + c/a + (1/a - c/a)*psi(c/a);

xx = 1e-2:0.01:10;
FuncVal = LLfunc(a,b,c,xx/ReducedArea);

% overlay histogram with expected fit
histogram(Areas(Areas<10),'Normalization','pdf','EdgeColor','none');
hold on
plot(xx,FuncVal/ReducedArea,'LineWidth',3);
hold off

% fit theoretical distribution
LLfunc = @(a,b,c,X) a*b^(c/a) / gamma(c/a) * X.^(c-1) .* exp(-b*X.^a);

KLScore=KLEstimator1DHybrid(log(Areas),UNIgamlog,1e-8)

% lets compare this with an analytical scoring mechanism on the
% histogrammed data
% KL-Divergence between two generalize gamma distributions
KLFunc = @(p1,p2) log(p1(1)*p2(2)^(-p2(3)/p2(1)) * gamma(p2(3)/p2(1)) / ...
    (p2(1)*p1(2)^(-p1(3)/p1(1))*gamma(p1(3)/p1(1)))) + ...
    (psi(p1(3)/p1(1))/p1(1) - 1/p1(1)*log(p1(2)) )*(p1(3)-p2(3)) ...
    + gamma( (p1(3)+p2(1))/p1(1)) / gamma( p1(3)/p1(1)) * (p1(2)^(-1/p1(1))/p2(2)^(-1/p2(1)))^(p2(1)) ...
    - p1(3)/p1(1);

arbCoeff = [1,1,1];
arbGamLog =  @(y) PDFgamlog(arbCoeff,y);

KLAnalytical = KLFunc(uniCoeff,arbCoeff)
KLvalidate = KLEstimator1DHybrid(log(Areas),arbGamLog,1e-8)






