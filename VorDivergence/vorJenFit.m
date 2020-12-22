% vorJenFit.m
% function to fit Jenny's data to the generalized gamma distribution and
% calculate the resulting KL-Divergence.
% Author: PKR, UPENN March 2019

% parse files from ui f
[files1, dataLoc1, indx1] = uigetfile('.csv', 'Select the csvs for treatment1', 'MultiSelect', 'on');
aa1=size(files1,2);

[files2, dataLoc2, indx2] = uigetfile('.csv', 'Select the csvs for treatment2', 'MultiSelect', 'on');
aa2=size(files2,2);

%options.WindowStyle = 'normal';
%str = inputdlg(...
%            'Name the output file.','Aggregate',[1 68],{''},options);
%str = str{1};

% parameters for the uniform distribution in 3D
a_uni = 1.16788;
b_uni = 4.04039;
c_uni = 4.79803;
uniCoef = [a_uni, b_uni, c_uni];

startGuess = [0.05,200,10];

spectral_val = 0.6;
%fitOptions = optimset('PlotFcns',@optimplotfval);
fitOptions.MaxIter = 1000;
fitOptions.TolFun = 1e-6;

%% likelihood function for the generalized gamma distribution
Lfunc = @(p,X) p(1)*p(2)^(p(3)/p(1)) / gamma(p(3)/p(1)) * X.^(p(3)-1) .* exp(-p(2)*X.^p(1));

%% negative log likelihood of the generalized gamma distribution
nLLfunc = @(p,X) -(p(3)/p(1))*log(p(2)) - log(p(1)) + gammaln(p(3)/p(1)) ...
    - (p(3)-1)*log(X) + p(2)*X.^p(1);

costFunc = @(p,X) -(p(3)/p(1))*log(p(2)) - log(p(1)) + gammaln(p(3)/p(1)) ...
    +mean(-(p(3)-1)*log(X) + p(2)*X.^p(1));


%% entropy function
EntropyVal = @(p) log(p(2)^(-1/p(1))*gamma(p(3)/p(1))/p(1)) + p(3)/p(1) + (1/p(1) - p(3)/p(1))*psi(p(3)/p(1));

%% KL Divergence from p1 to p2
KLFunc = @(p1,p2) log(p1(1)*p2(2)^(-p2(3)/p2(1)) * gamma(p2(3)/p2(1)) / ...
    (p2(1)*p1(2)^(-p1(3)/p1(1))*gamma(p1(3)/p1(1)))) + ...
    (psi(p1(3)/p1(1))/p1(1) - 1/p1(1)*log(p1(2)) )*(p1(3)-p2(3)) ...
    + gamma( (p1(3)+p2(1))/p1(1)) / gamma( p1(3)/p1(1)) * (p1(2)^(-1/p1(1))/p2(2)^(-1/p2(1)))^(p2(1)) ...
    - p1(3)/p1(1);

%% Experimental Data
%% Condition 1
% load the data sets, get coordinates per image
list1 = {};
xyz1 = {};
for iii=1:aa1
    
    data = fullfile(dataLoc1,files1{iii});
    fid = fopen(data);
    
    %% for Vutara formatted csv
    header1 = textscan(fid, '%s', 36, 'Delimiter',',');
    list1{iii} = textscan(fid, repmat('%f',1,36), 'HeaderLines',1,'Delimiter',',');
    xyz1{iii} = [list1{iii}{17}, list1{iii}{18}, list1{iii}{19}];
    fclose(fid);
    header1 = header1{1};
end

ExpVolumes1 = [];
VolVariances1 = zeros(aa1,1);
LocalVolumes1 = {};

for ii = 1:aa1
    % work on each file here
    x = xyz1{ii}(:,1);
    y = xyz1{ii}(:,2);
    z = xyz1{ii}(:,3);
    X = VoronoiVolumes([x,y,z],false);
    Volumes = X(:,4);
    % determine the perimeter of this file
    K = boundary(x,y,z,spectral_val);
    Keepers = true(length(Volumes),1);
    Keepers(K) = false;
    %sum(Keepers)/length(Keepers)
    ValidVolumes = Volumes(Keepers);
    % scale the areas so the expectation value is 1
    ValidVolumes = ValidVolumes/mean(ValidVolumes);
    LocalVolumes1{ii} = ValidVolumes;
    VolVariances1(ii) = var(ValidVolumes);
    ExpVolumes1 = [ExpVolumes1; ValidVolumes];
end

Expvolumes = [];
for ii = 1:aa1
    %if VolVariances1(ii) > 6
    %    continue
    %end
    ExpVolumes1 = [ExpVolumes1; LocalVolumes1{ii}];
end

% perform MLE estimation of Areas Distribution
%nLLsum1 = @(p) sum(nLLfunc(p,ExpVolumes1));
cost1 = @(p) costFunc(abs(p),ExpVolumes1);
%ContEst1 = fminsearch(nLLmean1,uniCoef,fitOptions);
%ContEst1 = fminsearch(cost1,uniCoef,fitOptions);
ContEst1 = fminsearch(cost1,startGuess,fitOptions);

%% Condition 2
% load the data sets, get coordinates per image
list2 = {};
xyz2 = {};
for iii=1:aa2
    
    data = fullfile(dataLoc2,files2{iii});
    fid = fopen(data);
    
    %% for Vutara formatted csv
    header2 = textscan(fid, '%s', 36, 'Delimiter',',');
    list2{iii} = textscan(fid, repmat('%f',1,36), 'HeaderLines',1,'Delimiter',',');
    xyz2{iii} = [list2{iii}{17}, list2{iii}{18}, list2{iii}{19}];
    fclose(fid);
    header2 = header2{1};
end

ExpVolumes2 = [];
VolVariances2 = zeros(aa2,1);
LocalVolumes2 = {};

for ii = 1:aa2
    % work on each file here
    x = xyz2{ii}(:,1);
    y = xyz2{ii}(:,2);
    z = xyz2{ii}(:,3);
    X = VoronoiVolumes([x,y,z],false);
    Volumes = X(:,4);
    % determine the perimeter of this file
    K = boundary(x,y,z,spectral_val);
    Keepers = true(length(Volumes),1);
    Keepers(K) = false;
    %sum(Keepers)/length(Keepers)
    ValidVolumes = Volumes(Keepers);
    % scale the areas so the expectation value is 1
    ValidVolumes = ValidVolumes/mean(ValidVolumes);
    LocalVolumes2{ii} = ValidVolumes;
    VolVariances2(ii) = var(ValidVolumes);
    ExpVolumes2 = [ExpVolumes2; ValidVolumes];
end

Expvolumes = [];
for ii = 1:aa2
    %if VolVariances2(ii) > 6
    %    continue
    %end
    ExpVolumes2 = [ExpVolumes2; LocalVolumes2{ii}];
end

% perform MLE estimation of Areas Distribution
%nLLmean2 = @(p) mean(nLLfunc(p,ExpVolumes2));
cost2 = @(p) costFunc(abs(p),ExpVolumes2);
%ContEst2 = fminsearch(nLLmean2,uniCoef,fitOptions);
%ContEst2 = fminsearch(cost2,uniCoef,fitOptions);
ContEst2 = fminsearch(cost2,startGuess,fitOptions);

%% Plotting Section
% going to hard-code bin size at 1000, remove volumes > 6
% get recommended binSize

%% hard-coded a maximum volume so that the histogram scales well
binNum = 1000;
figure;
hCon1 = histogram(ExpVolumes1(ExpVolumes1<6),binNum,'Normalization','pdf','EdgeAlpha',0);
hold on
hCon2 = histogram(ExpVolumes2(ExpVolumes2<6),binNum,'Normalization','pdf','EdgeAlpha',0);

% plot the fit coefficients of condition 1 to the histogram
%ContFuncVal1 = Lfunc(ContEst1,hCon1.BinEdges);
ContFuncVal1 = exp(-nLLfunc(ContEst1,hCon1.BinEdges));
plot(hCon1.BinEdges,ContFuncVal1,'c','LineWidth',3);
title('Histogrammed Volumes of Experimental Data and MLE fits','FontSize',14)
xlabel('Reduced Volumes','FontSize',14)
ylabel('Probability Density','FontSize',14)

% plot the fit coefficients of condition 2 to the histogram
%ContFuncVal2 = Lfunc(ContEst2,hCon2.BinEdges);
ContFuncVal2 = exp(-nLLfunc(ContEst2,hCon2.BinEdges));
plot(hCon2.BinEdges,ContFuncVal2,'r','LineWidth',3);
%title('Histogrammed Areas of Experimental Data and MLE fit','FontSize',14)
xlabel('Reduced Volumes','FontSize',14)
ylabel('Probability Density','FontSize',14)

legend('Condition 1','Condition 2','fit 1','fit 2')
hold off

%% Return the KL-Divergences

ContDiv1 = KLFunc(ContEst1,uniCoef);
ContDiv2 = KLFunc(ContEst2,uniCoef);

