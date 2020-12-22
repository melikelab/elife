% VorDivergence3D.m
% script to get KL-Divergence scores with uniform model in 3D

% get uniform distribution of KL-Divergence
% likelihood = log ab^(c/a)/Gamma(c/a) + (c-1)1/n sum^n_{i=1} log(x_i) -
% b/n sum^n_{i=1} log (x_i^a)

% simulation parameters
Nsamples = 1000;
Npoints = 400;
ReducedVolume = 2;
Dim = ReducedVolume*Npoints;

% statistical quantities
Volumes = [];
% generate points and extract areas over and over again
for ii = 1:Nsamples
    scale = Dim^(1/3);
    xyz = scale*rand(Npoints,3);
    % make the boundaries periodic
    xyz = [xyz; [xyz(:,1)+scale, xyz(:,2), xyz(:,3)]; ...
        [xyz(:,1)-scale, xyz(:,2), xyz(:,3)]; ...
        [xyz(:,1), xyz(:,2)+scale, xyz(:,3)]; ...
        [xyz(:,1), xyz(:,2)-scale, xyz(:,3)]; ...
        [xyz(:,1), xyz(:,2), xyz(:,3)+scale]; ...
        [xyz(:,1), xyz(:,2), xyz(:,3)-scale];];
    DT = delaunayTriangulation(xyz);
    [V,C] = voronoiDiagram(DT);
    VorDat = {V, C};
     
    % Calculate Volumes!
    X = zeros(Npoints,1);
    for pt=1:Npoints
        xt = V(C{pt},1);
        yt = V(C{pt},2);
        zt = V(C{pt},3);
        [~,X(pt)] = convhulln([xt yt zt],{'Qt','Qs'});
        %X(pt) = abs(sum( (xt([2:end 1]) - xt).*(yt([2:end 1]) + yt))*0.5);
    end
    Volumes = [Volumes; X(~isnan(X))];
end

% fit theoretical distribution
LLfunc = @(a,b,c,X) a*b^(c/a) / gamma(c/a) * X.^(c-1) .* exp(-b*X.^a);

% determine coefficients
a = 1.16788;
b = 4.04039;
c = 4.79803;

% Entropy of the generalized Gaussian distribution
EntropyVal = log(b^(-1/a)*gamma(c/a)/a) + c/a + (1/a - c/a)*psi(c/a);

xx = 1e-2:0.01:10;
FuncVal = LLfunc(a,b,c,xx/ReducedVolume);

% overlay histogram with expected fit
histogram(Volumes(Volumes<10),'Normalization','pdf','EdgeColor','none');
hold on
plot(xx,FuncVal/ReducedVolume,'LineWidth',3);
hold off

%% Voronoi plot
%figure;
%voronoi(xyz(:,1),xyz(:,2))
%[Vv,Cc] = voronoin(xyz, {'Qbb'});
%color = {'r' 'b' 'g' 'm' 'c' 'y' 'k' 'w'} ;
%for i = 1:length(Cc)
%   fill(Vv(Cc{i},1),Vv(Cc{i},2),char(randsample(color,1))) ;
%end

% define negative log likelihood function for losses
%nLLfunc = @(a,b,c,X) -(c/a)*log(b) - log(a) + gammaln(c/a) - (c-1)*log(X) + b*X.^a;

%Loss = nLLfunc(a,b,c,Volumes/ReducedVolume);

