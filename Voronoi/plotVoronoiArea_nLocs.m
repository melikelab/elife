

function plotVoronoiArea_nLocs(cluster,thresh,pix2nm,titlefilename)

%% initialize some params
mksz = 2;
% titlefilename = LL.filename( find(LL.filename==filesep,1,'last')+1:end-4);
titlefilename( titlefilename == '_' ) = '-';
nthresh = size(thresh,1);
colrs = lines(nthresh);
legstr = cell(nthresh,1);
h = nan(nthresh,1);

fg = figure;

sp = 1;
sb(sp) = subplot(1,1,1,'parent',fg);
hold(sb(sp),'on')
box(sb(sp),'on')

for t = 1:nthresh
    
    h(t) = plot(sb(sp), (pix2nm^2*cluster{t}.areas), cluster{t}.nLocs,...
        'Marker','s',...
        'LineStyle','none',...
        'Color',colrs(t,:),...
        'MarkerFaceColor',colrs(t,:),...
        'MarkerSize',mksz);
    
    legstr{t} = ['Thresh = ' num2str(thresh(t),'%.3g')];
    
end

ylabel(sb(sp),'Localizations per cluster')
xlabel(sb(sp),'Cluster Area (nm^{2})')
set(sb(sp),'XScale','log','YScale','log')
legend(h,legstr,'Location','NorthWest')
title(sb(sp),titlefilename)


% include file name as annotation
% titlefilename = filename;
% posA = [annotationXpos annotationYpos length(titlefilename(1:end-4))*0.02 0.03];
% if posA(3) > 1, posA(3) = 1;
% elseif posA(3) < 0, posA(3) = 0; end
% h = annotation(fg,'textbox',posA,...
%     'String',{titlefilename(1:end-4)},...
%     'FitBoxToText','on','EdgeColor','none',...
%     'FontSize',annotationFontSz);

end %of function
    