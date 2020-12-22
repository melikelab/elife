% vorShreyFit.m
% function to fit Shreyasi's data to the generalized gamma distribution and
% calculate the resulting KL-Divergence.

% list file locations
tree = 'O:\peter\ORP1L_mcherry_300nm';
control = 'control_.0156_16_300nm_data1.mat';
lovastatin = 'lovastatin_.0156_16_300nm_data1.mat';
u1866a = 'u1866a_.0156_16_300nm_data1.mat';

% parameters for the uniform distribution
a_uni = 1.07950;
b_uni = 3.03226;
c_uni = 3.31122;
uniCoef = [a_uni, b_uni, c_uni];

% function for the generalized gamma distribution
Lfunc = @(p,X) p(1)*p(2)^(p(3)/p(1)) / gamma(p(3)/p(1)) * X.^(p(3)-1) .* exp(-p(2)*X.^p(1));

% function for the negative log likelihood of the generalized gamma
% distribution
nLLfunc = @(p,X) -(p(3)/p(1))*log(p(2)) - log(p(1)) + gammaln(p(3)/p(1)) ...
    - (p(3)-1)*log(X) + p(2)*X.^p(1);

% entropy function
EntropyVal = @(p) log(p(2)^(-1/p(1))*gamma(p(3)/p(1))/p(1)) + p(3)/p(1) + (1/p(1) - p(3)/p(1))*psi(p(3)/p(1));

% load the three data sets
cont = load(fullfile(tree,control));
lova = load(fullfile(tree,lovastatin));
u18 = load(fullfile(tree,u1866a));

%% Control Data
ControlAreas = [];
ContLoopSize = size(cont.data1,1);

for ii = 1:ContLoopSize
    % work on the control data set here
    x = cont.data1{ii,1};
    y = cont.data1{ii,2};
    Area = cont.data1{ii,3};
    
    % determine the perimeter of this file
    K = boundary(x,y);
    Keepers = true(length(Area),1);
    Keepers(K) = false;
    ValidAreas = Area(Keepers);
    % scale the areas so the expectation value is 1
    ValidAreas = ValidAreas/mean(ValidAreas);
    ControlAreas = [ControlAreas; ValidAreas];
end

% perform MLE estimation of Areas Distribution
nLLsum = @(p) sum(nLLfunc(p,ControlAreas));
ContEst = fminsearch(nLLsum,uniCoef);

% get recommended binSize
[MbestCont,MAPcont] = binEstimator(ControlAreas);
figure;
hControl = histogram(ControlAreas,MbestCont,'Normalization','pdf','EdgeAlpha',0);
hold on

% plot the fit coefficients to the histogram
ContFuncVal = Lfunc(ContEst,hControl.BinEdges);
plot(hControl.BinEdges,ContFuncVal,'r','LineWidth',3);
title('Histogrammed Areas of Control Data and MLE fit','FontSize',14)
xlabel('Reduced Areas','FontSize',14)
ylabel('Probability Density','FontSize',14)
legend('histogram','fit distribution')
hold off

%% LovaStatin Data
StatinAreas = [];
StatinLoopSize = size(lova.data1,1);

for ii = 1:StatinLoopSize
    % work on the control data set here
    x = lova.data1{ii,1};
    y = lova.data1{ii,2};
    Area = lova.data1{ii,3};
    
    % determine the perimeter of this file
    K = boundary(x,y);
    Keepers = true(length(Area),1);
    Keepers(K) = false;
    ValidAreas = Area(Keepers);
    % scale the areas so the expectation value is 1
    ValidAreas = ValidAreas/mean(ValidAreas);
    StatinAreas = [StatinAreas; ValidAreas];
end

% perform MLE estimation of Areas Distribution
nLLsum = @(p) sum(nLLfunc(p,StatinAreas));
LovaEst = fminsearch(nLLsum,uniCoef);

% get recommended binSize
[MbestLova,MAPLova] = binEstimator(StatinAreas);
figure;
hLova = histogram(StatinAreas,MbestLova,'Normalization','pdf','EdgeAlpha',0);
hold on

% plot the fit coefficients to the histogram
LovaFuncVal = Lfunc(LovaEst,hLova.BinEdges);
plot(hLova.BinEdges,LovaFuncVal,'r','LineWidth',3);
title('Histogrammed Areas of LovaStatin Data and MLE fit','FontSize',14)
xlabel('Reduced Areas','FontSize',14)
ylabel('Probability Density','FontSize',14)
legend('histogram','fit distribution')
hold off

%% U1866a Data
U1Areas = [];
U1LoopSize = size(u18.data1,1);

for ii = 1:U1LoopSize
    % work on the control data set here
    x = u18.data1{ii,1};
    y = u18.data1{ii,2};
    Area = u18.data1{ii,3};
    
    % determine the perimeter of this file
    K = boundary(x,y);
    Keepers = true(length(Area),1);
    Keepers(K) = false;
    ValidAreas = Area(Keepers);
    % scale the areas so the expectation value is 1
    ValidAreas = ValidAreas/mean(ValidAreas);
    U1Areas = [U1Areas; ValidAreas];
end

% perform MLE estimation of Areas Distribution
nLLsum = @(p) sum(nLLfunc(p,U1Areas));
U1Est = fminsearch(nLLsum,uniCoef);

% get recommended binSize
[MbestU1,MAPU1] = binEstimator(U1Areas);
figure;
hU1 = histogram(U1Areas,MbestU1,'Normalization','pdf','EdgeAlpha',0);
hold on

% plot the fit coefficients to the histogram
U1FuncVal = Lfunc(U1Est,hU1.BinEdges);
plot(hU1.BinEdges,U1FuncVal,'r','LineWidth',3);
title('Histogrammed Areas of U1866a Data and MLE fit','FontSize',14)
xlabel('Reduced Areas','FontSize',14)
ylabel('Probability Density','FontSize',14)
legend('histogram','fit distribution')
hold off

% plot all fit curves on one plot
figure;
hold on
plot(hControl.BinEdges,ContFuncVal,'k','LineWidth',3);
plot(hLova.BinEdges,LovaFuncVal,'r','LineWidth',3);
plot(hU1.BinEdges,U1FuncVal,'c','LineWidth',3);
legend('Control','Lova Statin', 'U1866a');
hold off

% plot all step curves in one plot
figure;
hold on
[hC, Cedges] = histcounts(ControlAreas,MbestCont,'Normalization','pdf');
Cmids = (Cedges(1:end-1)+Cedges(2:end))/2;
plot(Cmids,hC,'k','LineWidth',3)

[hL, Ledges] = histcounts(StatinAreas,MbestLova,'Normalization','pdf');
Lmids = (Ledges(1:end-1)+Ledges(2:end))/2;
plot(Lmids,hL,'r','LineWidth',3)

[hU, Uedges] = histcounts(U1Areas,MbestU1,'Normalization','pdf');
Umids = (Uedges(1:end-1)+Uedges(2:end))/2;
plot(Umids,hU,'c','LineWidth',3)
legend('Control','Lova Statin', 'U1866a');
hold off

%% plot fits as CDFs
CDFit = @(p,x) gammainc(x.^(p(1))*p(2), p(3)/p(1));
Support = 0:0.01:6;
ContCDF = CDFit(ContEst,Support);
LovaCDF = CDFit(LovaEst,Support);
U1aCDF = CDFit(U1Est,Support);
figure;
hold on
plot(Support,ContCDF,'k','LineWidth',3);
plot(Support,LovaCDF,'r','LineWidth',3);
plot(Support,U1aCDF,'c','LineWidth',3);
legend('Control','Lova Statin', 'U1866a');
title('CDF of Fits')
ylabel('CDF')
xlabel('Reduced Areas')
ax2 = axes('Position',[0.4 0.2 0.4 0.4],'Box','on');
XX = 0:0.01:1;
CzCDF = CDFit(ContEst,XX);
LzCDF = CDFit(LovaEst,XX);
UzCDF = CDFit(U1Est,XX);
gca
plot(XX,CzCDF,'k','LineWidth',3)
hold on
plot(XX,LzCDF,'r','LineWidth',3)
plot(XX,UzCDF,'c','LineWidth',3)
hold off


%% Return the KL-Divergences
KLFunc = @(p1,p2) log(p1(1)*p2(2)^(-p2(3)/p2(1)) * gamma(p2(3)/p2(1)) / ...
    (p2(1)*p1(2)^(-p1(3)/p1(1))*gamma(p1(3)/p1(1)))) + ...
    (psi(p1(3)/p1(1))/p1(1) - 1/p1(1)*log(p1(2)) )*(p1(3)-p2(3)) ...
    + gamma( (p1(3)+p2(1))/p1(1)) / gamma( p1(3)/p1(1)) * (p1(2)^(-1/p1(1))/p2(2)^(-1/p2(1)))^(p2(1)) ...
    - p1(3)/p1(1);

ContDiv = KLFunc(ContEst,uniCoef)
LovaDiv = KLFunc(LovaEst,uniCoef)
U1Div = KLFunc(U1Est,uniCoef)

RevContDiv = KLFunc(uniCoef,ContEst)
RevLovaDiv = KLFunc(uniCoef,LovaEst)
RevU1Div = KLFunc(uniCoef,U1Est)

meanFunc = @(p) p(2)^(-1/p(1))*gamma((p(3)+1)/p(1))/(gamma(p(3)/p(1)));
modeFunc = @(p) p(2)^(-1/p(1))*( (p(3)-1)/p(1))^(1/p(1));
varFunc = @(p) p(2)^(-2/p(1)) * (gamma((p(3)+2)/p(1))/(gamma(p(3)/p(1))) ...
    - (gamma((p(3)+1)/p(1))/gamma(p(3)/p(1)))^2)

uniMean = meanFunc(uniCoef)
uniMode = modeFunc(uniCoef)
uniVar = varFunc(uniCoef)

contMean = meanFunc(ContEst)
empContMean = mean(ControlAreas)
contMode = modeFunc(ContEst)
[mxC,Cix]=max(hControl.Values);
empContMode = hControl.BinEdges(Cix)+hControl.BinWidth/2
contVar = varFunc(ContEst)
empContVar = var(ControlAreas)

lovaMean = meanFunc(LovaEst)
empLovaMean = mean(StatinAreas)
lovaMode = modeFunc(LovaEst)
[mxL,Lix]=max(hLova.Values);
empLovaMode = hLova.BinEdges(Lix)+hLova.BinWidth/2
lovaVar = varFunc(LovaEst)
empLovaVar = var(StatinAreas)

u1Mean = meanFunc(U1Est)
empU1Mean = mean(U1Areas)
u1Mode = modeFunc(U1Est)
[mxU,Uix]=max(hU1.Values);
empU1Mode = hU1.BinEdges(Uix)+hU1.BinWidth/2
u1Var = varFunc(U1Est)
empU1Var = var(U1Areas)


