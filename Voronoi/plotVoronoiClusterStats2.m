
function fig_h = plotVoronoiClusterStats2(numberLocs, clusterArea, nndXY, titlefilename)

titleFontSz = 9;
annotationFontSz = 10.5;
annotationYpos = 0.965;
annotationXpos = 0.38;
axisLabelFontSz = 8.5;
axisScaleFontSz = 8;
ylimIncrease = 1.05;
ylabelstr = 'Probability';
xaxLim = [1 1e3;...
          10 20e3;...
          0 500];

% fig_h = figure('OuterPosition',[285 0 600 900],'units','normalized');
fig_h = figure('Position',[150 150 1000 800]);
nr = 2;
nc = 2;
sp = 1;
%     nDomains = size(ClusterResults,1);

% include file name as annotation
% titlefilename = filename;
titlefilename( titlefilename == '_' ) = '-';
posA = [annotationXpos annotationYpos length(titlefilename(1:end-4))*0.02 0.03];
if posA(3) > 1, posA(3) = 1;
elseif posA(3) < 0, posA(3) = 0; end
h = annotation(fig_h,'textbox',posA,...
    'String',{titlefilename(1:end-4)},...
    'FitBoxToText','on','EdgeColor','none',...
    'FontSize',annotationFontSz);

% distribution of number of localizations per nanodomain
sb(sp) = subplot(nr,nc,sp,'Parent',fig_h,'XScale','log','XMinorTick','on');
box on
plotdata = numberLocs;
% Jan 4 2016 update for plotting with stair function
% [bins,ct] = plotStairs(plotdata);
% stairs(sb(sp),bins,ct,...
%     'Color','b',...
%     'LineStyle','-');
% xlabel(sb(sp),'Number of Localizations/nanodomain','FontSize', axisLabelFontSz)
% ylabel(sb(sp),ylabelstr,'FontSize', axisLabelFontSz)
% title(sb(sp),['Median # Loc/nanodomain = ' num2str( median(plotdata),'%.3g' ) ],...
%     'FontSize', titleFontSz)
% axis(sb(sp),[ xaxLim(sp,:) 0 max(ct)*ylimIncrease])
% set(sb(sp),'XScale','log','XLim',[1 2000])
% set(sb(sp),'XScale','log','XTickLabel',{'1' '10' '100' '1000'},'FontSize',axisScaleFontSz)
% % adjust annotation position
% posP = get(sb(1),'Position');
% posA = get(h,'Position');
% posA(1) = posP(1)+0.5*(posP(3)-posA(3));
% if posA(1) < 0, posA(1) = 0; end
% posA(2) = annotationYpos;
% set(h,'Position',posA)
% sp = sp+1;
% 
% % distribution of nanodomain area
% % areas = [pi*results(:,4).*results(:,5), pi*results(:,6).^2, pi*results(:,7).^2].*original_pixel_size^2;
% % ellipse area, circle radius=avg(sigX,sigY), circle radius = sqrt(sigX^2+sigY^2)
% sb(sp) = subplot(nr,nc,sp,'Parent',fig_h,'XScale','log','XMinorTick','on');
% box on
% plotdata = clusterArea;
% % Jan 4 2016 update for plotting with stair function
% [bins,ct] = plotStairs(plotdata);
% stairs(sb(sp),bins,ct,...
%     'Color','b',...
%     'LineStyle','-');
% xlabel(sb(sp),'Nanodomain Area (nm^2)','FontSize', axisLabelFontSz)
% ylabel(sb(sp),ylabelstr,'FontSize', axisLabelFontSz)
% title(sb(sp),['Nanodomain Area: median = ' num2str( median(plotdata),'%.3g' ) ...
%     ', average = ' num2str( mean(plotdata),'%.3g' ) ],...
%     'FontSize', titleFontSz)
% axis(sb(sp),[ xaxLim(sp,:) 0 max(ct)*ylimIncrease])
% set(sb(sp),'XScale','log','XTickLabel',{'10' '100' '1000' '10000'},'FontSize',axisScaleFontSz)
% sp = sp+1;
% 
% % distribution of nearest neighbor distances
% sb(sp) = subplot(nr,nc,sp,'Parent',fig_h,'XScale','linear','XMinorTick','on');
% box on
% plotdata = nndXY;
% % Jan 4 2016 update for plotting with stair function
% [bins,ct] = plotStairs(plotdata,[],'linear');
% stairs(sb(sp),bins,ct,...
%     'Color','b',...
%     'LineStyle','-');
% xlabel(sb(sp),'Nearest Neighbor Distance (nm)','FontSize', axisLabelFontSz)
% ylabel(sb(sp),ylabelstr,'FontSize', axisLabelFontSz)
% title(sb(sp),['Nearest Neighbor Distance: median = ' num2str( median(plotdata),'%.3g' ) ...
%     ', average = ' num2str( mean(plotdata),'%.3g' ) ' nm'],...
%     'FontSize', titleFontSz)
% axis(sb(sp),[ xaxLim(sp,:) 0 max(ct)*ylimIncrease])
% sp = sp+1;
% 
% % Relation between custer area and number of localizations
% sb(sp) = subplot(nr,nc,sp,'Parent',fig_h,...
%     'XScale','log',...
%     'XMinorTick','on',...
%     'YScale','log',...
%     'YMinorTick','on'...
%     );
% box on
% colr = lines(1);
% plot(sb(sp), clusterArea, numberLocs ,...
%     'Marker','s',...
%     'LineStyle','none',...
%     'Color',colr,...
%     'MarkerFaceColor',colr,...
%     'MarkerSize',4)
% ylabel(sb(sp),'Localizations per cluster')
% xlabel(sb(sp),'Cluster area (nm^2)')
% set(sb(sp),'XScale','log','YScale','log')
% end % of function