%plotDensityScript.m
% script to generate density data for voronoi analysis

[fn, path]=uigetfile('*.bin','Select the Voronoi Analysed .bin files(s)','MultiSelect','on');
if ischar(fn)
    fn = {fn};
end
nfiles=size(fn,2);

tree = path;

for i = 1:nfiles
data = fn{i};

% load the data
LL = Insight3( fullfile( tree, data ) );
xy = LL.getXYcorr;  % get xy positions in pixels

%pix2nm = 116;
%barWidth = 1;
%fig_h = plotVoronoiDiag_mapArea(xy,  pix2nm, barWidth);

%inputs
%1: coordinates of your data
%2: lower and upper threshold of the area value in the coordinates input
%3: conversion factor of native coordinates to nm
%4: not needed, but if you have the voronoi data, you can put it there
%5: the size of the scale bar
fig_h = plotVoronoiDensities(xy,[2.5e-3,0.04],98.69,[],1000);
fig_h.Color='black';
set(gcf,'InvertHardCopy','off');

savefile = fullfile(path, [fn{i}(1:end-15) '.png']);
saveas(gcf,savefile,'png')
end