function [fig_h, h] = plotVoronoiMCdat(Histograms, thresh, signif)
%% plot Voronoi area graph output from VoronoiMonteCarlo_JO
fig_h = figure('Name','Voronoi Clustering Threshold');
hold on
box on
mxfactor = 1.05;
countCols = [2,size(Histograms,2)];
ylimTop = ceil(mxfactor*max(...
        [Histograms(:,countCols(1));Histograms(:,countCols(2))]...
        ));
% base = [2 5]; % for plotting the filled area without clipping the outline box
base = [min(Histograms(:,1))*0.05 min(Histograms(:,2))*0.05];
clear area

if ~isempty(thresh)
    mfact = 2.25;
    colr = lines(8)*mfact;
    colr = colr([6 5 2 4 7],:);
    colr(colr>1) = 1;
    colr = repmat(colr,ceil(length(thresh)/size(colr,1)),1);
    harea = nan(length(thresh),1);
    legstrA = cell(1,length(thresh));
    for c = 1:length(thresh)
    harea(c) = ...
        area([base(1),thresh(c)],[1,1]*ylimTop-base(2),base(2),...intersection(2),...
        'LineStyle','none',...
        'FaceColor',colr(c,:));
    legstrA{c} = ['Thresh = ' num2str(thresh(c),'%.2g')];
    end
end
colr = {[0 0.4 0.75];... blue
        [0 0.7 0];... green
        [0.8 0.3 0]};% orange
colr{4} = colr{3};
h = nan(4,1);
if size(Histograms,2)==5 %legacy from October 2016
    for c = 2:5
        h(c-1) = plot(Histograms(:,1), Histograms(:,c),...
            'Color',colr{c-1},...
            'LineWidth',1.4);
    end
    xmax = max(Histograms(:,1));
else % updated January 25, 2017 since MonteCarlo_x values ~= data_x
    c = 2;% data column
    h(c-1) = plot(Histograms(:,1), Histograms(:,c),...
        'Color',colr{c-1},...
        'LineWidth',1.4);
    for c = 4:6
        h(c-2) = plot(Histograms(:,3),Histograms(:,c),...
            'Color',colr{c-2},...
            'LineWidth',1.4);
    end
    xmax = max([Histograms(:,1);Histograms(:,3)]);
end
% axis([0 Inf, 0 ylimTop])
set(gca,'XLim',[0 xmax],'YLim',[0 ylimTop])
xlabel('Voronoï polygon area')
ylabel('Count')
title('Distribution of Voronoï polygon area')
legstr = {'Data','Monte Carlo',[num2str(signif) '% conf bounds']};
if ~isempty(thresh)
    legstr = [legstr,legstrA];
    h = [h(1:3);harea];
end
legend(h,legstr,...
    'Location','NorthEast')

end % of function
