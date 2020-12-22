% J.J.Otterstrom, Matlab 2013b
% Program designed to receive an Insight3 .bin molecule list acquried
% during a z-calibration measurement and generate the z-calibration fit
% results for Wx and Wy.  Results are both plotted and saved
%
% INPUT
% acq = acquisition parameters of the z-calibration images
%     [scan_range_nm, step_size_nm]
%     If not input, the user will be asked to define it
% ML = path & name of molecule list file, i.e.: path\filename.bin
%     Optional input, otherwise user must select it
% Optional Variable Input Pairs (varargin)
% 'stageMovement' = Optional string input describing stage movement during
%     data acquisition
%     values are: 'ascending' (default, i.e. -600 to +600 nm) 
%                 'decending' (i.e. +600 to -600 nm) BUT this is
%                     automatically updated when using calibration data
%                     acquired on STORM3, since the actual stage
%                     positions are used.
% 'mSigma' = multiplication factor 'm' for determining where the
%     threshold of m*sigma that the residual of the linear approximation 
%     to the fit line must fall below; sigma is the standard deviation of
%     the residuals between the fit line and the mean valued ratio Wx/Wy
%     data used to generate the fit line.

% OUTPUT
% structure having the following fields:
%     .pathname = string of the molecule list path location
%     .filename = string of the molecule list file name
%     .xyzwidth = [x_coord, y_coord, z_coord, Wx_val, Wy_val, Wx/Wy_val]
%     .meanvals = [z, mean_Wx, mean_Wy, mean_ratio_Wx/Wy], these are the
%         values used to obtain fit parameters
%     .fitparam = [W_ox, z_rx, g_x, B_x, A_x;
%                  W_oy, z_ry, g_y, B_y, A_y]
%     .fit = [z, fit_Wx, fit_Wy, fit_ratio_Wx/Wy]


function Output = fitZcalib(acq,MLname,varargin)

% UPDATES
% - Feb 10, 2015 - changed Wx and Wy definitions
% - Feb 23, 2015 - reverted to original Wx and Wy definitions
%     > Proper use of Insight3 with the appropriate microscope ensures
%     the calibration curve is correct
%         >> OriginalSTORM microscope = I3 from 2011.05.26
%         >> STORM3 microscopes = I3 from 2014.11.06 (Joe's version)
% - May 21, 2015 - fixed error during the definition of variable 'Fr'
%     arising when there was a frame without detections
% - May 28, 2015 - updating to include 'stageMovement' optional input and
%     to automatically extract real stage position information stored in
%     the .inf files acquired on the STORM3 microscope. This involved
%     updating the DAX.m function by Joe.
% - June 19, 2015 - Included functionality to estimate the linear regime
%     of the Wx/Wy ratio fit line. This is performed by performing a
%     taylor expansion around each z-depth between -450 <= z <= +450 nm 
%     Then the z-range where the residual from the taylor expansion is
%     about the same as that for the calibration fit to the mean data 
%     values is determined and set as the "linear regime"
% - December 2015/March 2016 - included the 'stageMovement' variable to 
%     define the movement of the stage as "ascending" or "descending". When
%     the objective moves, then the input value should be reversed, since
%     their movements are opposed.
% - May 18, 2016 - changed y-axis of plots to fixed limits
% - July 14, 2016 - updated subfunction processargs() to avoid variable
%     definition errors in static workspaces. 
%     Also changed default stage movement to 'decending' since STORM3 with
%     the piezo-stage follows this movement


%% unpack variable input arguments and define defaults when appropriate
optargin = size(varargin,2);
if rem(optargin,2)~=0
    error('invalid syntax. Input property name followed by value');
end

defaultVals = { 'decending';   ... stageMovement
                 1.5;          ... mSigma
              };
[stageMove,mSigma] ...
                = processargs(varargin,optargin,defaultVals);
            
%% Testing values
% stageMovement.value = 'ascending';
% stageMovement.userInput = false

%% check for correct input scan and step sizes
if nargin == 0 || ~isvector(acq) || max(size(acq)) ~= 2
    accept = 0;
    while accept == 0
        disp('Z-calibration acquisition parameters not established')
        acq(1) = input('Z-scan range in nanometers (nm) ==> ');
        acq(2) = input('Scan step size in nanometers (nm) ==> ');
        accept = input(['Input values: [' num2str(acq(1),'%.4g') ', ' ...
            num2str(acq(2),'%.2g') '], is this correct? 0=no, 1=yes  ==> ']);
    end
end    

%% Load Molecule list

tol = 1e-6; % tolerance for accepting complex valued output from Wfit fitting
war = 0; % index for counting the warnings that may appear 

% check for input molecule list name or ask user to locate it
if sum(nargin == [0,1])==1 || isempty(strfind(MLname,'\')) || sum(strcmp(MLname(end-3:end),{'.bin'}))==0
    % query user to select .bin molecule list file saved by Insight3
    [filename, pathname] = uigetfile({'*.bin'},'Select Insight3 Molecule List');
    ML = Insight3([pathname filename]);
else
    if exist(MLname,'file')
        ML = Insight3(MLname);
        idx = strfind(MLname,'\');
        pathname = MLname(1:idx(end));
        filename = MLname(idx(end)+1:end);
    else
        error(['No file found with name: ' MLname])
    end
end


%% Begin algorithm

% define the z-calibration curve to be fit: 
% Zcalib = @(p,x)( Wo*(sqrt(1+((x-g)/z)^2+a*((x-g)/z)^3+b*((x-g)/z)^4)) );
% p(1)=Wo  p(2)=z  p(3)=g  p(4)=b   p(5)=a
Wfit = @(p,x)( p(1).*sqrt(1+...
                              ((x-p(3))./p(2)).^2+...
                        p(5).*((x-p(3))./p(2)).^3+...
                        p(4).*((x-p(3))./p(2)).^4 ...
                           ) );

% remove file type from filename variable
filename = filename(1:end-4);

% get frame numbers and corresponding z-values
nFr = ML.numFrames;
% obtain z-values corresponding to each frame
%  actual stage position values are available from STORM3 .inf files and
%  accessible through the DAX.m function
% 1) parse filename to find what the .dax file associated with the molecule
%    list, presuming that the root of the ML file name is the same
w=0;
while exist([ML.filename(1:end-w) '.dax'],'file') ~= 2
    w = w+1;
    if w >= length(ML.filename)
        w = nan(1,1);
        break
    end
end
if ~isnan(w)
    dax = DAX([ML.filename(1:end-w) '.dax']);
    % 2) here extract the actual stage position data
    [stageData,header] = dax.getZCalStagePos();
    % 3) check to see if the data is there, i.e. STORM3 calibration curve
else
    stageData = [];
    header = [];
end
if ~isempty(stageData)
    if ~isempty(strfind(header{3},'um'))
%     if abs(stageData(1,3)<100) % stage positions in um
        multfc = 1000; % convert um to nm
%     else
    elseif ~isempty(strfind(header{3},'nm'))
        multfc = 1; % stage position is in nm
    end
    steps = multfc*(stageData(1,3) - stageData(:,3));
    % check to make sure that steps is ascending (default)
    tstfr = 2;
    while tstfr < size(stageData,1)
        if steps(tstfr) ~= 0 
            if steps(tstfr) > 0 && ~stageMove.userInput % do not adjust value if the user provided the input
                stageMove.value = 'ascending'; % else stick with default/input
            end
            break
        end
    	tstfr = tstfr+1;
    end
else
    % not STORM3, so approximate the z-position as traditionally done
    steps = transpose(0:acq(2):acq(1));
end

%% define the direction of stage movement parameter
switch stageMove.value
    case 'ascending'
        movement = 1;
    case 'decending'
        movement = -1;
end

if nFr ~= size(steps,1) % NSTORM calibration condition
    % Need to remove first and last 20 frames (40 in total)
    ibgn = find(ML.data(:,14)==21,1,'first');
    iend = find(ML.data(:,14)==nFr-20,1,'last');
    % re-define the molecule list to only include these values
    ML.data = ML.data(ibgn:iend,:);
    nLoc = iend-ibgn+1;
    nFr = nFr-40;
    
else % Labview calibration condition, no extraneous frames    
    nLoc = ML.getNumberOfMolecules();
end

% Fr = [unique(ML.data(:,14)), steps]; % [frame#, z-value]
Fr = [unique(ML.data(:,14)), movement*steps(unique(ML.data(:,14)))];
nFr = size(Fr,1);

% Convert the eliptical gaussian fit parameters into Wx and Wy values
% [Wx,Wy] = calcWidths(ML.data);
[Wx,Wy] = calcWidths(ML);
ratio = [Wx./Wy, ML.data(:,14)]; % [Wx/Wy_ratio, frame number]

% determine which frame corresponds to z = 0
mnvals=zeros(nFr,3); % [mean_Wx, mean_Wy, Mean_ratio
for i = 1:nFr    
    idx = find(ratio(:,2) ==  Fr(i,1));
    mnvals(i,:) = [mean(Wx(idx)) mean(Wy(idx)) mean(ratio(idx,1))];     
end

% [~,z0idx] = min(abs(mnvals(:,3)-1));
[~,z0idx(1)] = min(abs(mnvals(:,3)-1)); 
[~,z0idx(2)] = min(abs(mnvals(:,1)-mnvals(:,2)));
% expect that the best index to use for z=0 is closest to the middle of the
% recording
[~,i]=min(Fr(z0idx,1)-nFr/2);
% Fr = [Fr, Fr(:,2)-Fr(z0idx,2)];
z0idx = z0idx(i);
Fr = [Fr, Fr(:,2)-Fr(z0idx,2)];

if max(Fr(:,3)) < 450 
    war = war+1;
    Output.warning{war} = 'Elongation does not extend above +450 nm, may give poor fits';
    warning(Output.warning{war})
end
if min(Fr(:,3)) > -450
    war = war+1;
    Output.warning{war} = 'Elongation does not extend below -450 nm, may give poor fits';
    warning(Output.warning{war})
end  
%% set each Wx and Wy to a z value, then fit & evaluate the linear region
% [Xc, Yc, Z, Wx, Wy, Wx/Wy]
xyzwidths = zeros(nLoc,6);
for i = 1:nFr
    
    idx = find(ratio(:,2) ==  Fr(i,1));
    xyzwidths(idx,:) = [ML.data(idx,[3 4]), Fr(i,3)*ones(size(idx,1),1), Wx(idx,:), Wy(idx,:), ratio(idx,1)];
    
end

% set the initial guess of the fit parameters
% Wo is ~300, z is ~300, a is ~±0.1, b is ~±0.01
pinitX = [300, 300, NaN, 0.1, 0.01]; pinitY = pinitX;
% if Wx/Wy < 1 when z < 0, then gx < 0 and gy > 0
idx1 = find(xyzwidths(:,3) == Fr(1,3));
idx2 = find(xyzwidths(:,3) == Fr(2,3));
ratioInit = mean(xyzwidths([idx1;idx2],6));
if ratioInit < 1
    pinitX(3) = -400*movement;
    pinitY(3) = 400*movement;
else
    pinitX(3) = 400*movement;
    pinitY(3) = -400*movement;
end

% Perform non-linear fitting
resid = zeros(nLoc,2); % W_x & W_y residuals
[FitParam(1,:),resid(:,1)] = nlinfit(xyzwidths(:,3),xyzwidths(:,4),Wfit,pinitX); % Wx
[FitParam(2,:),resid(:,2)] = nlinfit(xyzwidths(:,3),xyzwidths(:,5),Wfit,pinitY); % Wy

% Define output data
Output.pathname = pathname;
Output.filename = filename;
Output.xyzwidth = xyzwidths;
Output.meanvals = [Fr(:,3), mnvals];

% check to see if fitting results are real valued
if ~isempty( find(imag(FitParam)~=0,1) ) % then some parameters have complex values
    % check to see if the imaginary contribution is above a tolerance
    if ~isempty( find( abs(imag(FitParam)./real(FitParam)) >= tol ,1) )
        error('Complex Valued Fit Parameter Beyond Accepted Tolerance: skip data set, increase tolerance or figure out why it''s complex-valued')
    else
        if ~isempty( find(imag(FitParam(1,:))~=0,1) )
            war = war+1;
            Output.warning{war} = 'Wx fit parameters were returned with complex values';
            warning(Output.warning{war})
        end
        if ~isempty( find(imag(FitParam(2,:))~=0,1) )
            war = war+1;
            Output.warning{war} = 'Wy fit parameters were returned with complex values';
            warning(Output.warning{war})
        end
        Output.fitparamOrig = FitParam;
        % redefine parameters using the real part, since the imaginary part is
        % negligible
        FitParam = real(FitParam);
    end
end
    
% prepare fit parameters and the resulting calibration fit for output
Output.fitparam = FitParam;
Output.fit(:,1) = Fr(:,3);
for i = 1:3
    if i ~= 3
        Output.fit(:,i+1) = Wfit(real(FitParam(i,:)),Output.fit(:,1));
    else
        Output.fit(:,i+1) = Output.fit(:,2)./Output.fit(:,3);
    end
end

%% find a region where linear fitting of Wx/Wy is appropriate

% Taylor Expansion to estimate the linear region of the fit line
% I) define derivative values needed for taylor expansion & linar approximation
fitCalib = @(fp,z)(Wfit(fp(1,:),z)./Wfit(fp(2,:),z));
dWfit_dz = @(p,x)( (p(1)/(2*p(2))) .* ...
                    (... 
                        1+((x-p(3))./p(2)).^2+p(5).*((x-p(3))./p(2)).^3+p(4).*((x-p(3))./p(2)).^4 ...
                    ).^(-0.5) .* ...                    ) .* ...
                    (...
                        2.*((x-p(3))./p(2)) + ...
                        3.*p(5).*((x-p(3))./p(2)).^2 + ...
                        4.*p(4).*((x-p(3))./p(2)).^3 ...
                    ) );
dfitCalib_dz = @(fp,z)( ...
                -1.*Wfit(fp(1,:),z).*(Wfit(fp(2,:),z)).^(-2).*dWfit_dz(fp(2,:),z) + ...
                (Wfit(fp(2,:),z)).^(-1).*dWfit_dz(fp(1,:),z) ...
                );

% II) Make linear Taylor approx & find where closest to fit & mean-value data
zrangeIdx = intersect( find( Output.meanvals(:,1) > -450), ...
                       find( Output.meanvals(:,1) < 450) );

% nSigma = 1.5; % 2*sigma = 95% of Normal distrib, 3*sigma = 99% of normal distrib
sigma = std(Output.meanvals(:,4) - Output.fit(:,4),1);
sigmaThresh = mSigma*sigma;
% d= cell(length(zrangeIdx),1);
LinFit.all = nan(length(zrangeIdx),6);
clear idx
% format for output
Output.LinearRange.DataColVal = {'z_low', 'z_high', 'expansion-center', 'slope', ...
                                 'intercept', 'z-range', 'R^2', 'p-Value', ' Rank_score'};
for i = 1:length(zrangeIdx)
    % point where taylor expansion is performed
    a = Output.meanvals(zrangeIdx(i),1);
    slope = dfitCalib_dz(FitParam,a);
    intercept = fitCalib(FitParam,a);
    % calculate fit line's linear taylor expansion approximation & residuals
    taylorExp = intercept+(Output.fit(zrangeIdx,1)-a).*slope;
    residTaylor = abs(Output.fit(zrangeIdx,4) - taylorExp);
    % find the z-range where the residual from the taylor expansion is
    % about the same as that for the calibration fit to the mean data values
    closeIdx = find(residTaylor <= sigmaThresh);
    % find the number of continuous z-depths where the linear, taylor
    % expansion is a good approximation of the fit.
    d = diff(closeIdx);
%     d{i} = diff(closeIdx);
    idx(1) = find(d==1,1,'first');
    idx(2) = find(d==1,1,'last'); 
    tstidx = find(d(idx(1):idx(2)) > 1)+idx(1)-1;
    %
    while ~isempty(tstidx)
        % move either the idx(1) up or idx(2) down
        [~,ii] = min(abs(idx-tstidx(1)));
        if ii == 1
            idx(1) = tstidx(1)+1;
        elseif ii == 2
            idx(2) = tstidx(1)-1;
        end
        tstidx = find(d(idx(1):idx(2)) > 1)+idx(1)-1;
    end
    idx(2) = idx(2)+1; % since the diff() operation is d(1:end-1)=closeIdx(2:end)-closeIdx(1:end-1)
    % Save info to output
    LinFit.all(i,[1 2]) = Output.fit(zrangeIdx(closeIdx(idx)),1);
    if LinFit.all(i,1) > LinFit.all(i,2)
        LinFit.all(i,[1 2]) = Output.fit(zrangeIdx(closeIdx(idx([2 1]))),1);
    end
    
    LinFit.all(i,3) = a; % center of taylor expansion
    LinFit.all(i,4) = slope;
    LinFit.all(i,5) = intercept;
    % calculate z-range    
    LinFit.all(i,6) = LinFit.all(i,2) - LinFit.all(i,1); % z-range
    % calculate R^2 
    r2 = R2_calc( ...
                 Output.fit(zrangeIdx(closeIdx(idx(1):idx(2))),4), ...
                 taylorExp(closeIdx(idx(1):idx(2))) ...
         );
    LinFit.all(i,7) = r2;
    % calculate the likeness of the fit & linear approx to the mean data
    data = Output.meanvals(zrangeIdx(closeIdx(idx(1):idx(2))),4);
    residFit = data - Output.fit(zrangeIdx(closeIdx(idx(1):idx(2))),4);
    residLin = data - taylorExp(closeIdx(idx(1):idx(2)));
    [~,pVal] = kstest2(residFit,residLin);
    LinFit.all(i,8) = pVal;
    LinFit.all(i,9) = pVal*r2;

end

Output.LinearRange.allData = LinFit.all;

% 3) find where pVal > 0.05 such that there is no statistical difference 
%    between the fit line & linear approx. If none, then take the biggest 
%    pVal and include a warning that it may be a crap calib
LinFit.good = LinFit.all( LinFit.all(:,8)>0.05, :); 
if isempty(LinFit.good)
    linZrange = BestLinFit( LinFit.all );
    war = war+1;
    Output.warning{war} = ['Maxmimum p-value for calib fit & linear similiarity = ' ...
            num2str(linZrange(8),'%.3g') '; calibration untrustworthy'];
    warning(Output.warning{war})
    
elseif size(LinFit.good,1) > 1 % more than one 'best'
    Output.LinearRange.DataGood = LinFit.good;
    linZrange = BestLinFit( LinFit.good );
else
    linZrange = LinFit.good;
end

Output.LinearRange.best = linZrange;

%% Plot: data + mean values + fit curves
nrow = 3;
ncol = 3;

fontsize = 10;

pink = [1 0.72 0.72];
ltbl = [0.72 0.81 1];
gray = ones(1,3)*0.8;
ltcolr = [pink;ltbl;gray];

rd = [0.85 0 0];
bl = [0 0.43 0.93];
gy = [0.5 0.5 0.5];
dot = [rd;bl;gy];

% for linear fit
yel = [0.95 0.55 0];

z = xyzwidths(:,3); % z-position for each localization
zz = Fr(:,3);  % z-position corresponding to each frame
% Output.fit(:,1) = Fr(:,3); 
yy = zeros(size(zz,1),2);

plotOrder = {1;4;[2,3];[5,6]};
ylims1 = [200,1200;...
         200,1200;...
         0.45 2.5];
dataOrder = [4,5,6];
titleData = {'W_x','W_y','W_x/W_y','W_x & W_y'};
% yup = [1300;1300;2.5];
zlims = [-450 450];

fh = figure('OuterPosition',[500 150 900 900]);
orient(fh,'landscape')


% plot Wx, Wx and Wx/Wy each vs. Z
for i = 1:3
    
    ydata = xyzwidths(:,dataOrder(i));
    c = ltcolr(i,:);
    
    subplot(nrow,ncol,plotOrder{i},'Parent',fh)
    hold on
    box on
    % plot raw data lightly in background
    plot(z,ydata,'s','Color',c,'MarkerFaceColor',c,'MarkerSize',4)
    % plot average values at each z-value
    plot(zz,mnvals(:,i),'o','Color',dot(i,:),'MarkerFaceColor',dot(i,:),'MarkerSize',7)
    yy(:,i) = Output.fit(:,i+1);
    h = nan(2,1);
    % plot fit line
    h(2) = plot(zz,yy(:,i),'k-','LineWidth',2.5);
    
    if i ~= 3
        xlabel('Z-position (nm)','FontSize',fontsize)
        axis([-600 600 200 1200])
    end
    ylabel(titleData{i},'FontSize',fontsize)
    if i == 3
        % plot best linear fit
        h(1) = plot(linZrange(1:2),(linZrange(1:2)-linZrange(3))*linZrange(4)+linZrange(5),...
            '--','Color',yel,'LineWidth',4.5);
        % set axis limits
%         [~,i3] = min(abs(z-zlims(1))); [~,i4] = min(abs(z-zlims(2)));
        axis([zlims, ...
            ylims1(i,:)])%max([0 min(ydata(min([i3 i4]):max([i3 i4])))]) min([yup(i) max(ydata(min([i3 i4]):max([i3 i4])))]) ])
%         axis([ zlims, 0.9*max([0 min(mnvals(i1:i2,3))]), 1.15*min([yup(i) max(mnvals(i1:i2,3))]) ])
        filename(filename == '_') = '-';
        title({[titleData{i} ' vs. z: ' filename];...
            ['Linear Z-range estimate: ' num2str(linZrange(1),'%.3g') ...
                'nm to ' num2str(linZrange(2),'%.3g') 'nm, linear R^2 = ' ...
                num2str(linZrange(7),'%.3g')]} ...
            ,'FontSize',fontsize)
        legend(h,'Linear Regime','W_x/W_y Fit','Location','North')
    else
        switch stageMove.value
            case 'ascending'
                axis([zz(1)-10 zz(end)+10 ylims1(i,:)])% max([0 min(ydata)]) min([yup(i) max(ydata)])])
            case 'decending'
                axis([zz(end)-10 zz(1)+10 ylims1(i,:)])% max([0 min(ydata)]) min([yup(i) max(ydata)])])
        end
        title([titleData{i} ' vs. z'],'FontSize',fontsize)
    end

end

%% Plot Wx & Wy vs. Z, also save the calibration fit parameters to a text
% file and display it for easy copy/paste on the command screen
j = 4;
subplot(nrow,ncol,plotOrder{j},'Parent',fh)
hold on
box on
legstr = cell(2,1); subtitle = cell(2,1);
axordr = {'x','y'};
h = zeros(4,1);
m = 0;
% print file name to command screen
fprintf(['\nFit Params for ' filename '\n'])
% open file for writing the calibration fit parameters
fid = fopen([pathname filename '_ZcalibFitParam.txt'],'w+');
for i = 1:2
    
    p = FitParam(i,:);
    m = m+1;
    % plot mean values
    h(m)=plot(zz,mnvals(:,i),'o','Color',dot(i,:),'MarkerFaceColor',dot(i,:),'MarkerSize',7);
    yy(:,i) = Wfit(p,zz);
    m = m+1;
    h(m) = plot(zz,yy(:,i),'k-','LineWidth',2.5);
    % plot fit lines
%     axis([ zlims, 0.9*max([0 ylim(1)]), 1.15*min([yup(i) ylim(2)]) ])
    if i == 2
        axis([ zlims, 0, 1200]);%, 1.07*max([mnvals(:,1);mnvals(:,2)]) ]);%, 1.15*min([yup(i) ylim(2)]) ])
    end
    legstr{i} = titleData{i};
    % Create parameter output into subplot title
    subtitle{i} = ['W_{' axordr{i} 'o} = ' num2str(p(1),'%.4g') ...
                 ', z_' axordr{i} ' = ' num2str(p(2),'%.4g') ...
                 ', g_' axordr{i} ' = ' num2str(p(3),'%.4g') ...
                 ', B_' axordr{i} ' = ' num2str(p(4),'%.4g') ...
                 ', A_' axordr{i} ' = ' num2str(p(5),'%.4g')];
    % format fit parameters for command line display & writing to file
    Wstr = ['w' axordr{i} '0=' sprintf('%6.2f',p(1)) ';'];
    zstr = ['zr' axordr{i} '=' sprintf('%6.2f',p(2)) ';'];
    gstr = ['g' axordr{i} '=' sprintf('%6.2f',p(3)) ';'];
    if p(4) < 0, Bstr = ['B' axordr{i} '=' sprintf('%7.5f',p(4)) ';'];
    else Bstr = ['B' axordr{i} '=' sprintf('%6.5f',p(4)) ';']; end
    if p(5) < 0, Astr = ['A' axordr{i} '=' sprintf('%7.5f',p(5)) ';'];
    else Astr = ['A' axordr{i} '=' sprintf('%6.5f',p(5)) ';']; end
    fprintf([Wstr zstr gstr '\nC' axordr{i} '=0.0000;' Bstr Astr '\n']) 
    fprintf(fid,[Wstr zstr gstr '\nC' axordr{i} '=0.0000;' Bstr Astr '\n']);
end
fprintf('\n')
fclose(fid);
order = [1 3];
xlabel('Z-position (nm)','FontSize',fontsize)
ylabel(titleData{j},'FontSize',fontsize)
% title({subtitle{1}; subtitle{2} },'FontSize',10)
title([subtitle{1} '; ' subtitle{2}],'FontSize',7)
legend(h(order),legstr,'Location','North')

%% Plot symetricity measurements
% First reverse Wy fit and plot on top of Wx fit
sp1 = subplot(nrow,ncol,7,'Parent',fh);
hold(sp1,'on')
box(sp1,'on')
fp = [1 -1];
lnstl = {'-.','-'};
legadd = {'',' flipped'};
ylims2 = [200 1000];
clear idx
hp = zeros(2,1);
for i = 1:2 
    hp(i) = plot(sp1,fp(i)*Output.meanvals(:,1),Output.meanvals(:,i+1),...
        'o','Color',dot(i,:),'MarkerFaceColor',dot(i,:),'MarkerSize',7);
end
for i = 1:2
    ydata = Output.fit(:,i+1);
    [~,idx(1)] = min(abs(Output.fit(:,1)+425));
    [~,idx(2)] = min(abs(Output.fit(:,1)-425));
    if idx(2) < idx(1), idx = sort(idx); end
    plot(sp1,fp(i)*Output.fit(:,1),ydata,'Color',ltcolr(i,:),...
        'LineStyle',lnstl{i},'LineWidth',2.5)
    ylims2 = [ min([ylims2(1) min( ydata(idx(1):idx(2)) )]) max([ylims2(2) max( ydata(idx(1):idx(2)) )]) ];
    legstr{i} = [titleData{i} legadd{i}];
end
axis(sp1,[-400 400 ylims2])
xlabel(sp1,'Z-position (nm)','FontSize',fontsize)
ylabel(sp1,'Width Value (nm)','FontSize',fontsize)
hl = legend(hp,legstr,'Location','North');
set(hl,'FontSize',fontsize-1)
title(sp1,'W_x & W_y symmetricity')

% Second Plot the inverse of the Wx/Wy ratio when that ratio is below 1
% LinReg = 200; % half-range where ratio vs. Z is approximated to be linear
sp2 = subplot(nrow,ncol,[8,9],'Parent',fh);
hold(sp2,'on')
box(sp2,'on')
% find where the ratio is above or below one
idxBig = cell(3,1); idxSml = idxBig;
idxBig{1} = find(Output.meanvals(:,4) > 1);
idxBig{2} = find(Output.fit(:,4) > 1);
idxSml{1} = find(Output.meanvals(:,4) <= 1);
idxSml{2} = find(Output.fit(:,4) <= 1);
% indicies for setting the y-axis limits
clear idx
[~,idx(1)] = min(abs(Output.fit(:,1)+300));
[~,idx(2)] = min(abs(Output.fit(:,1)-300));
if idx(2) < idx(1), idx = sort(idx); end
% plot ratio >= 1
plot(sp2,Output.meanvals(idxBig{1},1),Output.meanvals(idxBig{1},4),...
    'o','Color',dot(3,:),'MarkerFaceColor',dot(3,:),'MarkerSize',7)
plot(sp2,Output.fit(idxBig{2},1),Output.fit(idxBig{2},4),...
    '-','Color','k','LineWidth',2.5)
% plot ratio < 1 as 1/ratio
plot(sp2,Output.meanvals(idxSml{1},1),1./Output.meanvals(idxSml{1},4),...
    'o','Color',dot(3,:),'MarkerFaceColor',dot(3,:),'MarkerSize',7)
plot(sp2,Output.fit(idxSml{2},1),1./Output.fit(idxSml{2},4),...
    '-','Color','k','LineWidth',2.5)
% plot the linear fit
x = linZrange(1):acq(2):linZrange(2);
y = (x-linZrange(3))*linZrange(4)+linZrange(5);
plot(sp2,x, abs(y-1)+1,...
    '--','Color',yel,'LineWidth',3.5)
% x1 = Output.fit(idxSml{2}(Output.fit(idxSml{2},1)>=linZrange(1)),1);
% plot(sp2,x1, abs( linZrange(4)*(x1-linZrange(3))+linZrange(5)-1 )+1,...
%     '--','Color',yel,'LineWidth',3.5)
% x2 = Output.fit(idxBig{2}(Output.fit(idxBig{2},1)<=linZrange(2)),1);
% plot(sp2,x2, abs( linZrange(4)*(x2-linZrange(3))+linZrange(5)-1 )+1,...
%     '--','Color',yel,'LineWidth',3.5)
% set axes labels & ranges
xlabel(sp2,'Z-position (nm)','FontSize',fontsize)
ylabel(sp2,'Abs( (W_x/W_y-1) )','FontSize',fontsize)
axis(sp2,[-400 400 0.95 2.4])% 1.1*max([ Output.fit(idx,4);1./Output.fit(idx,4)])])
set(sp2,'YTick',1:0.2:2.4)

%% subfunction to define variable input variables
    function [stageMovementInput,mSigma] = ...
            processargs(argsin,noptargin,defaults)
        
        stageMovement = [];
        mSigma = [];
        
        for idxnoptarg = 1:2:noptargin
            if ~ischar(argsin{idxnoptarg})
                error('invalid syntax. Input property name followed by value');
            end
            eval([argsin{idxnoptarg} ' = varargin{idxnoptarg+1};'])
        end
        
%         % initialize variable
%         stageMovement = struct();
%         stageMovement.value = defaults{1};
%         stageMovement.userInput = false;
        
        % check for the stage movement definition
        if ~exist('stageMovement','var') || isempty(stageMovement)
            stageMovementInput.value = defaults{1};
            stageMovementInput.userInput = false;
        else
            vals = strfind({'ascending','decending'},stageMovement);
            if isempty(vals{1}) && isempty(vals{2})
                in = [];
                while isempty(in)
                    in = input('Invalid input for Stage Movement. Type ''1'' for ascending or ''2'' for decending ==> ');
                    if in == 1
                        stageMovementInput.value = 'ascending';
                        stageMovementInput.userInput = true;
                    elseif in == 2
                        stageMovementInput.value = 'decending';
                        stageMovementInput.userInput = true;
                    else
                        in = [];
                    end
                end
            else
                stageMovementInput.value = stageMovement;
                stageMovementInput.userInput = true;
            end
        end
        
        % check for the mSigma multiplicative factor determining the
        % linear-fit residual threshold
        if ~exist('mSigma','var') || isempty(mSigma)
            mSigma = defaults{2};
        end
        
    end % end of processargs

% %% subfunction to calculate the widths of each fit
%     function [Wx,Wy] = calcWidths(ML)
%         
%         % Input: 18-column Insight3 molecule list, as read by 'readInsight3.m'
%         % Output: The widths in 'x' and 'y' for each fitted localization
%         
%         W = ML(:,7);
%         phi = ML(:,8);
%         ax = ML(:,9);
%         
%         Ws = W./sqrt(ax);
%         Wl = W.*sqrt(ax);
%         
%         % use this definition if phi is in radians
%         % Wx = Ws.*cos(phi);
%         % Wy = Wl.*cos(phi);
%         % I found out that phi is in degrees
%         % The following Wx,Wy definition used until Feb 10, 2015
%         % It was mistakenly altered until Feb 23, 2015 to the definitions
%         % below: Wx(Wl) and Wy(Ws).  
%         % Instead, the definitions using Wx(Ws) and Wy(Wl) are correct.
%         % The update was reversed because the confusion arose from proper
%         % use of the correct version of Insight 3
%         % OriginalSTORM microscope must use the 2011.05.26 version of I3
%         % NSTORM and STORM3 microscopes must use the 2014.11.06 version of I3
%         Wx = Ws.*cosd(phi);
%         Wy = Wl.*cosd(phi);
%         % Updated on Feb 10, 2015 because I found out that the definition
%         % was reversed
%         % These two definitions removed Feb 23, 2015
% %         Wy = Ws.*cosd(phi);
% %         Wx = Wl.*cosd(phi);
%         
%     end % end of calcWidths sub-function

%% subfunction to calculate the R2 values
    function [r2, r2adj] = R2_calc(obs,model,df)
        
        % df = degrees of freedom for the input model
        if ~exist('df','var'), df = 1; end
        n = length(obs);
        df_t = n-1;
        df_e = n-df-1;
        
        SStot = sum( (obs - mean(obs)).^2 );
        
        SSerr = sum( (obs - model).^2 );
        
        r2 = 1 - SSerr/SStot;
        
        if r2 < 0, r2 = 0; end
        
        if nargout > 1
            r2adj = 1 - (SSerr/df_e)/(SStot/df_t);
            if r2adj < 0, r2adj = 0; end
        end
    end % end of R2_calc sub-function

%% subfunction to extract the best linear fit for asymmetric z-ranges
    function best = BestLinFit( LinFits )
% Output.LinearRange.DataColVal = {'z_low', 'z_high', 'expansion-center', 'slope', ...
%                                  'intercept', 'z-range', 'R^2', 'p-Value', ' Rank_score'};
        % remove all entries that do not contain z=zero as this is
        % necessarily the most sensitive region
        for jj = 1:size(LinFits,1)
            if LinFits(jj,1) < 0 && LinFits(jj,2) > 0 
                % do nothing, all is normal
            else % does not contain zero
                % reset the rank score to zero to prevent it being chosen
                LinFits(jj,9) = 0;
            end
        end
        % take the one that has the highest rank factor: p-value & r^2-value
        [~,Lidx] = max( LinFits(:,9) );
        if length(Lidx) > 1
            % take the largest z-range if there are multiple 'max' p-values
            [~,idxx] = max( LinFits(Lidx,6) );
            if length(idxx) > 1
                % take the most symmetric option available
                [~,idxxx] = abs( mean( LinFits(Lidx(idxx),[1 2]) ,2) );
                if length(idxxx) > 1
                    % just take the first one, since it appears not to matter
                    Lidx = Lidx(idxx(idxxx(1)));
                else
                    Lidx = Lidx(idxx(idxxx));
                end
            else
                Lidx = Lidx(idxx);
            end
        end
        best = LinFits(Lidx,:);
        
    end % end of BestLinFit sub-function


end
% end of function