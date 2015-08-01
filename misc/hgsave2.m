function hgsave(varargin)
% HGSAVE  Saves an HG object hierarchy to a MAT file.
%
% HGSAVE(H, 'Filename') saves the objects identified by handle array H
% to a file named 'Filename'.  If 'Filename' contains no extension,
% then the extension '.fig' is added.  If H is a vector, none of the
% handles in H may be ancestors or descendents of any other handles in
% H.
%
% HGSAVE('Filename') saves the current figure to a file named
% 'Filename'.
%
% HGSAVE(..., 'all') overrides the default behavior of excluding
% non-serializable objects from the file.  Items such as the default
% toolbars and default menus are marked as non-serializable, and are
% normally not saved in FIG-files, because they are loaded from
% separate files at figure creation time.  This allows the size of
% FIG-files to be reduced, and allows for revisioning of the default
% menus and toolbars without affecting existing FIG-files. Passing
% 'all' to HGSAVE insures non-serializable objects are also saved.
% Note: the default behavior of HGLOAD is to ignore non-serializable
% objects in the file at load time, and that can be overridden using
% the 'all' flag with HGLOAD.
% HGSAVE(..., 'all') produces a warning and the option will be removed in a
% future release.
%
% HGSAVE(..., '-v6') saves a FIG-file that can be loaded by versions
% prior to MATLAB 7. When creating a figure to be saved and used in a
% version prior to MATLAB 7 use the 'v6' option to the plotting
% commands. See the help for PLOT, BAR and other plotting commands for
% more information.
%
% See also HGLOAD, SAVE.

%   Copyright 1984-2009 The MathWorks, Inc.
%   D. Foti  11/10/97

error(nargchk(1, inf, nargin,'struct'));

% Get the handle to save and the file to save to.
[h, filename, varargin] = localGetHandleAndFile(varargin);

% Process the remaining arguments
[SaveAll, SaveOldFig, SaveFlags] = localParseOptions(varargin);

% Default fig format values
hgS = struct(...
    'type', {}, ...
    'handle', {}, ...
    'properties', {}, ...
    'children', {}, ...
    'special', {});
hgO = [];
        
% Decide which save code path to use
if 0
    % Warn if user passed in 'all' flag
    if SaveAll
        warning( 'MATLAB:hgsave:DeprecatedOption', ...
            'The ''all'' option to hgsave will be removed in a future release.');
    end
    
    hgS = hgsaveStructDbl(h, SaveAll);
    SaveVer = '070000';
    SaveOldFig = true;
else
    % Warn if user passed in 'all' flag
    if SaveAll
        warning( 'MATLAB:hgsave:DeprecatedOption', ...
            'The ''all'' option to hgsave has been removed.');
    end
    
    if SaveOldFig
        hgS = hgsaveStructClass(h);
        SaveVer = '070000';
    else
        hgO = hgsaveObject(h);
        SaveVer = '070000';
    end
end

% Revision encoded as 2 digits for major revision,
% 2 digits for minor revision, and 2 digits for
% patch revision.  This is the minimum revision
% required to fully support the file format.
% e.g. 070000 means 7.0.0
SaveVars.(sprintf('hgS_%s', SaveVer)) = hgS;
if ~SaveOldFig
    SaveVars.(sprintf('hgO_%s', SaveVer)) = hgO;
end

save(filename, '-struct', 'SaveVars', SaveFlags{:});



function [h, filename, args] = localGetHandleAndFile(args)
% Look for a handle and/or filename in input arguments

% pull off handle + 'filename,' or just 'filename'
if ischar(args{1})
    h = gcf;
    filename = args{1};
    args(1) = [];
else
    h = args{1};
    if any(~ishghandle(h))
        E = MException('MATLAB:hgsave:InvalidHandle','Invalid handle');
        E.throwAsCaller();
    end
    filename = args{2};
    args(1:2) = [];
end

% Add an implicit ".fig" to the filename if no extension is specified
[path, file, ext] = fileparts(filename);

% fileparts returns everything from the last . to the end of the
% string as the extension so the following test will catch
% an extension with 0, 1, or infinity dots.
% for example, all these filenames will have .fig added to the end:
%  foo.
%  foo..
%  foo.bar.
%  foo.bar...

if isempty(ext) || strcmp(ext, '.')
    filename = fullfile(path, [file , ext, '.fig']);
end


function [SaveAll, SaveOldFig, SaveFlags] = localParseOptions(args)
% Parse optional flags from remaining arguments

if ~iscellstr(args)
    E = MException('MATLAB:hgsave:InvalidOption','Options must all be strings.');
    E.throwAsCaller();
end

% Default values
SaveAll = false;
SaveOldFig = true;   % Default is old-style fig files for now.
SaveFlags = {};

for n = 1:length(args)
    opt = args{n};
    if strcmpi(opt, 'all')
        SaveAll = true;
    elseif strcmpi(opt, '-figv7')
        SaveOldFig = true;
    elseif ~isempty(regexp(opt, '^-v[\d.]+$', 'once', 'start'))
        % -vX.Y flag is passed on to the save function
        SaveFlags = {opt};
    else
        % Error on any unrecognised option
        E = MException('MATLAB:hgsave:InvalidOption','Invalid option: %s.', opt);
        E.throwAsCaller();
    end
end

function hgS = hgsaveStructClass(h)
%hgsaveStructClass Save object handles to a structure.
%
%  hgsaveStructClass converts handles into a structure ready for saving.
%  This function is called when MATLAB is using objects as HG handles.

%   Copyright 2009 The MathWorks, Inc.

hgS = handle2struct(h);
