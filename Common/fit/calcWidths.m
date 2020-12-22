% J.Otterstrom, Matlab 2013b
% For input list of localizations determined from Insight3 and read into
% MATLAB using Joe Borbely's Insight3.m class, this function will return
% the width in X and Y (Wx and Wy) for a 3D list of localizations.
%
% function call:
% [Wx,Wy] = calcWidths(ML)
%   units of Wx and Wy are in nanometers
function [Wx,Wy] = calcWidths(ML)


if ~isa(ML,'Insight3')
    error('Input must be an Insight3 localization list ... champ')
    
else
    
    % Input: 18-column Insight3 molecule list, as read by 'readInsight3.m'
    % Output: The widths in 'x' and 'y' for each fitted localization
    W = ML.getColumn('width');
    phi = ML.getColumn('phi');
    ax = ML.getColumn('aspect');
%     
%     W = ML(:,widthcol);
%     phi = ML(:,phicol);
%     ax = ML(:,aspectcol);
%     
    Ws = W./sqrt(ax);
    Wl = W.*sqrt(ax);
    
    % use this definition if phi is in radians
    % Wx = Ws.*cos(phi);
    % Wy = Wl.*cos(phi);
    
    Wx = Ws.*cosd(phi);
    Wy = Wl.*cosd(phi);
    
end
end