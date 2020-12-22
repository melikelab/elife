% function files = AskForFiles(extensions, start_folder, window_title)
%
% Display a dialog window for selecting 1 or more files.
%
% Inputs
% ------
% extensions (optional) : string or cell/char array
%   The files that match these extensions are visible in the dialog window
%   Examples
%           'dax'
%           '.dax'
%           '*.dax'
%           'sample*.dax'
%           ['dax';'bin']
%           {'dax';'dcimg'}
%           {'sam*.dax';'H2_*.dcimg';'bin'}
%   If not set, then the default is 'All Files (*.*)'
%
% start_folder (optional) : string  
%   The folder that the dialog window start in
%   If not set, then the default is the Current Folder of Matlab
%
% window_title (optional) : string
%   The text to display in the titlebar of the dialog window
%   If not set, then the default title is 'Select file(s)'
%
% Returns
% -------
% A cell of file names. If no files were selected then the cell is empty.
%
function files = AskForFiles(extensions, start_folder, window_title)

default_extensions = '*.*';
default_folder = '';
default_title = 'Select file(s)';

switch nargin
    case 3
        FilterSpec = parseExtension(extensions);
        DefaultName = start_folder;
        DialogTitle = window_title;
    case 2
        FilterSpec = parseExtension(extensions);
        DefaultName = start_folder;
        DialogTitle = default_title;
    case 1
        FilterSpec = parseExtension(extensions);
        DefaultName = default_folder;
        DialogTitle = default_title;
    otherwise
        FilterSpec = default_extensions;
        DefaultName = default_folder;
        DialogTitle = default_title;
end

if exist(DefaultName, 'dir') ~= 7
    DefaultName = default_folder;
end

if ~isa(DialogTitle, 'char')
    DialogTitle = default_title;
end

[FileName, PathName, ~] = uigetfile(FilterSpec, DialogTitle, DefaultName, 'MultiSelect', 'on');

if isa(FileName, 'double')
    % then the user clicked the "Cancel" button on the UI
    files = {};
elseif isa(FileName, 'char')
    % then the user selected 1 file
    files = {fullfile(PathName, FileName)};
else
    % then the user selected multiple files
    files = cell(length(FileName), 1);
    for i=1:length(FileName)
        temp = fullfile(PathName, FileName(i));
        files{i} = temp{1};
    end
end

end

function out = parseExtension(input)
original = input;

if isa(input, 'char')
    % handle the case where the input is, for example, ['*.dax'; '*.bin']
    nrows = size(input, 1);
    temp = cell(nrows,1);
    for i = 1:nrows
        temp{i} = input(i,:);
    end
    input = temp;
end

assert(isa(input, 'cell'), 'The extensions parameter must be a string or cell')

out = cell(length(input),2);
for i=1:length(input)
    if isa(input{i}, 'char')
        asterix_index = strfind(input{i}, '*');
        period_index = strfind(input{i}, '.');
        if isempty(asterix_index) && isempty(period_index)
            out{i,1} = sprintf('*.%s', input{i});
            out{i,2} = upper(input{i});
        elseif isempty(asterix_index) && ~isempty(period_index)
            temp = strrep(input{i}, '.', '');
            out{i,1} = sprintf('*.%s', temp);
            out{i,2} = upper(temp);
        elseif length(asterix_index) == 1 && strcmp(input{i}(1), '*')
            temp = strrep(strrep(input{i}, '*', ''), '.', '');
            out{i,1} = sprintf('*.%s', temp);
            out{i,2} = upper(temp);
        elseif ~isempty(period_index)
            [~, name, ex] = fileparts(input{i});                    
            out{i,1} = [name ex];
            out{i,2} = upper(ex(2:end));
        else
            out = original;
            break;
        end        
    end
end

end

