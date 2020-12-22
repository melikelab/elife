function Output = fitZcalib(acq,MLname)

%%% INPUT
%%% acq = acquisition parameters of the z-calibration images
%%%     [scan_range_nm, step_size_nm]
%%%     If not input, the user will be asked to define it
%%% ML = path & name of molecule list file, i.e.: path\filename.bin
%%%     Optional input, otherwise user must select it

%%% OUTPUT
%%% structure having two fields:
%%%     .xyzwidth = [x_coord, y_coord, z_coord, Wx_val, Wy_val, Wx/Wy_val]
%%%     .fitparam = [W_ox, z_rx, g_x, B_x, A_x;
%%%                  W_oy, z_ry, g_y, B_y, A_y]

%%% UPDATES
%%% - Feb 10, 2015 - changed Wx and Wy definitions
%%% - Feb 23, 2015 - reverted to original Wx and Wy definitions
%%%     > Proper use of Insight3 with the appropriate microscope ensures
%%%     the calibration curve is correct
%%%         >> OriginalSTORM microscope = I3 from 2011.05.26
%%%         >> STORM3 microscopes = I3 from 2014.11.06 (Joe's version)

% % Testing
% %%
% pathname = 'V:\Jason_O\Data\141016 - 3D phantom\ECM+ beads 3D mult FOV\z-calibrations\';
% filename = 'zcalib_DualView_000_listGREEN.bin';
% acq = [1200, 25];
% %%
% ML = Insight3(['V:\Jason_O\LiveSTORM\141125 - WView test3\zcalib 141125\'...
%     'zcalib_1600_001_listGREEN.bin']);
% acq = [1600, 25];
% %%

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

% Testing:
% MLname = 'V:\Lara\neurons\141210 - zcalib\short_001_list.bin'; acq = [1200,25];

if sum(nargin == [0,1])==1 || isempty(strfind(MLname,'\')) || sum(strcmp(MLname(end-3:end),{'.bin'}))==0
    % query user to select .bin molecule list file saved by Insight3
    [filename, pathname] = uigetfile({'*.bin'},'Select Insight3 Molecule List');
    ML = Insight3([pathname filename]);
else
    ML = Insight3(MLname);
    idx = strfind(MLname,'\');
    pathname = MLname(1:idx(end));
    filename = MLname(idx(end)+1:end);
end

%% Begin algorithm

% define the z-calibration curve to be fit: 
% Zcalib = @(p,x)( Wo*(sqrt(1+((x-g)/z)^2+a*((x-g)/z)^3+b*((x-g)/z)^4)) );
% p(1)=Wo  p(2)=z  p(3)=g  p(4)=b   p(5)=a
Zcalib = @(p,x)( p(1).*sqrt(1+...
                              ((x-p(3))./p(2)).^2+...
                        p(5).*((x-p(3))./p(2)).^3+...
                        p(4).*((x-p(3))./p(2)).^4 ...
                           ) );

% remove file type from filename variable
filename = filename(1:end-4);

% get frame numbers and corresponding z-values
nFr = ML.numFrames;
steps = transpose(0:acq(2):acq(1));
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
Fr = [unique(ML.data(:,14)), steps]; % [frame#, z-value]

% Convert the eliptical gaussian fit parameters into Wx and Wy values
[Wx,Wy] = calcWidths(ML.data);
ratio = [Wx./Wy, ML.data(:,14)]; % [Wx/Wy_ratio, frame number]

% determine which frame corresponds to z = 0
mnvals=zeros(nFr,3); % [mean_Wx, mean_Wy, Mean_ratio
for i = 1:nFr    
    idx = find(ratio(:,2) ==  Fr(i,1));
    mnvals(i,:) = [mean(Wx(idx)) mean(Wy(idx)) mean(ratio(idx,1))];     
end

[~,z0idx] = min(abs(mnvals(:,3)-1));
Fr = [Fr, Fr(:,2)-Fr(z0idx,2)];

% set each Wx and Wy to a z value
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
    pinitX(3) = -400;
    pinitY(3) = 400;
else
    pinitX(3) = 400;
    pinitY(3) = -400;
end
% Perform non-linear fitting
resid = zeros(nLoc,2); % W_x & W_y residuals
[FitParam(1,:),resid(:,1)] = nlinfit(xyzwidths(:,3),xyzwidths(:,4),Zcalib,pinitX); % Wx
[FitParam(2,:),resid(:,2)] = nlinfit(xyzwidths(:,3),xyzwidths(:,5),Zcalib,pinitY); % Wy

% Define output data
Output.xyzwidth = xyzwidths;
Output.fitparam = FitParam;
Output.meanvals = [Fr(:,3), mnvals];

%% Plot: data + mean values + fit curves
nrow = 4;
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

z = xyzwidths(:,3);
zz = Fr(:,3); 
Output.fit(:,1) = Fr(:,3); 
yy = zeros(size(zz,1),2);

plotOrder = {1;4;[2,3];[5,6]};
dataOrder = [4,5,6];
titleData = {'W_x','W_y','W_x/W_y','W_x & W_y'};
yup = [1300;1300;3];
zlims = [-450 450];
i1 = find(zz==zlims(1),1,'first'); i2 = find(zz==zlims(2),1,'last');
ylim = [min([zlims(1);mnvals(i1:i2,1);mnvals(i1:i2,2)]) ...
        max([zlims(2);mnvals(i1:i2,1);mnvals(i1:i2,2)])];

fh = figure;
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
    if i ~= 3
        yy(:,i) = Zcalib(FitParam(i,:),zz);
    else
        yy(:,i) = yy(:,1)./yy(:,2);
    end
    Output.fit(:,i+1) = yy(:,i);
%     axis([zz(1)-10 zz(end)+10 max([0 min(ydata)]) min([yup(i) max(ydata)])])
    % plot fit line
    plot(zz,yy(:,i),'k-','LineWidth',2.5)
    
    if i ~= 3
        xlabel('Z-position (nm)','FontSize',fontsize)
    end
    ylabel(titleData{i},'FontSize',fontsize)
    if i == 3
        i3 = find(z==zlims(1),1,'first'); i4 = find(z==zlims(2),1,'last');
        axis([zlims, max( [0 min(ydata(i3:i4))] ) min( [yup(i) max(ydata(i3:i4))] ) ])
%         axis([ zlims, 0.9*max([0 min(mnvals(i1:i2,3))]), 1.15*min([yup(i) max(mnvals(i1:i2,3))]) ])
        filename(filename == '_') = '-';
        title([titleData{i} ' vs. z: ' filename],'FontSize',fontsize)
    else
        axis([zz(1)-10 zz(end)+10 max([0 min(ydata)]) min([yup(i) max(ydata)])])
        title([titleData{i} ' vs. z'],'FontSize',fontsize)
    end

end

% Plot Wx & Wy vs. Z, also save the calibration fit parameters to a text
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
    yy(:,i) = Zcalib(p,zz);
    m = m+1;
    h(m) = plot(zz,yy(:,i),'k-','LineWidth',2.5);
    % plot fit lines
    axis([ zlims, 0.9*max([0 ylim(1)]), 1.15*min([yup(i) ylim(2)]) ])
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

% Plot Residuals of Wx and Wy
zlims2 = [zz(1)-acq(2) zz(end)+acq(2)*1.5];
sp = subplot(nrow,ncol,[7,8,9],'Parent',fh);
hold on
box on
for i = 1:2
boxplot(sp,resid(:,i),z,...
    'plotstyle','compact',...    'medianstyle','line',...
    'colors',dot(i,:),...
    'outliersize',4,...
    'positions',zz+acq(2)*0.5*(i-1));
    if i == 1, set(gca,'XTickLabel',{' '}), end
end
plot(zlims2,[0 0],'k-','LineWidth',0.75)
axis([zlims2 [-1 1]*4*std([resid(:,1);resid(:,2)],1)])
xlabel('Z-position (nm)','FontSize',fontsize)
ylabel('Residual','FontSize',fontsize)
% legend('W_x','W_y')
title('W_x & W_y residuals','FontSize',fontsize)

% Plot Residuals of Wx/Wy
sp = subplot(nrow,ncol,[10,11,12],'Parent',fh);
hold on
box on
residratio = zeros(size(z,1),1);
for i = 1:size(zz,1)
    idx = find(z==zz(i));
    residratio(idx,1) = xyzwidths(idx,6)-yy(i,3);
end
boxplot(sp,residratio,z,...
    'plotstyle','compact',...    
    'colors',dot(3,:),...
    'outliersize',4,...
    'positions',zz)
plot(zlims2,[0 0],'k-','LineWidth',0.75)
axis([zlims2 [-1 1]*4*std(residratio,1)])
xlabel('Z-position (nm)','FontSize',fontsize)
ylabel('Residual','FontSize',fontsize)
title('W_x/W_y residuals','FontSize',fontsize)

%% subfunction to calculate the widths of each fit
    function [Wx,Wy] = calcWidths(ML)
        
        %%% Input: 18-column Insight3 molecule list, as read by 'readInsight3.m'
        %%% Output: The widths in 'x' and 'y' for each fitted localization
        
        W = ML(:,7);
        phi = ML(:,8);
        a = ML(:,9);
        
        Ws = W./sqrt(a);
        Wl = W.*sqrt(a);
        
        % use this definition if phi is in radians
        % Wx = Ws.*cos(phi);
        % Wy = Wl.*cos(phi);
        % I found out that phi is in degrees
        % The following Wx,Wy definition used until Feb 10, 2015
        % It was mistakenly altered until Feb 23, 2015 to the definitions
        % below: Wx(Wl) and Wy(Ws).  
        % Instead, the definitions using Wx(Ws) and Wy(Wl) are correct.
        % The update was reversed because the confusion arose from proper
        % use of the correct version of Insight 3
        % OriginalSTORM microscope must use the 2011.05.26 version of I3
        % NSTORM and STORM3 microscopes must use the 2014.11.06 version of I3
        Wx = Ws.*cosd(phi);
        Wy = Wl.*cosd(phi);
        % Updated on Feb 10, 2015 because I found out that the definition
        % was reversed
        % These two definitions removed Feb 23, 2015
%         Wy = Ws.*cosd(phi);
%         Wx = Wl.*cosd(phi);
        
    end
end