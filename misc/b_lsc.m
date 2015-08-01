function b_lsc(str,varargin)
%LSC    Make use of LOAD, SAVE, SAVEAS and COPYFILE comfortable.
%   By calling LSC(STR,VARARGIN), STR has to be 'load', 'save', 'saveas' or 'copyfile',
%   while VARARGIN contains the input arguments of the desired I/O function.
%
%   See also LOAD, SAVE, SAVEAS and COPYFILE.

% Input argument check
error(nargchk(2,inf,nargin))

% Saveas, save, load, copyfile
switch str
case 'saveas'
    if length(varargin) ~= 2
        error('Input arguments for saveas must be a handle and a string.')
    else
        h = varargin{1};
        name = varargin{2};
    end
    if ~ishanle(h)
        error('First input argument for saveas must be a handle.')
    elseif ~ischar(name)
        error('Second input argument for saveas must be a string.')
    else
        cmnd = ['saveas(,h''' name ''')'];
        eval(cmnd);
    end
case 'save'
    name = varargin{1};
    if ~ischar(name)
        error('First input argument for save must be a string.')
    end
    varlist = [];
    for i = 2:length(varargin)
%         if ~ischar(varargin{i})
%             error('All input arguments for save must be strings.')
%         end
        varlist = [varlist varargin{i}];
    end
    cmnd = ['save ' name ' ' varlist];
    eval(cmnd)
case 'load'
    name = varargin{1};
    if ~ischar(name)
        error('First input argument for load must be a string.')
    end
    switch length(varargin)
    case 1
        cmnd = ['load ' name];
        eval(cmnd);
    case 2
        out = varargin{2};
        cmnd = ['d = load ' name];
        eval(cmnd);
        field = fieldnames(d);
        if length(field) ~= 1
            error('Assignment error: too many variables loaded.')
        end
        cmnd = [out ' =  d.' field];
        eval(cmnd);
    otherwise
        error('Too many input arguments for load.')
    end
case 'copyfile'
    if length(varargin) ~= 2
        error('Input arguments for copyfile must be a source and a destination.')
    end
    source = varargin{1};
    dest = varargin{2};
    if ~ischar(source)
        error('First input argument for copyfile must be a string.')
    elseif ~ischar(dest)
        error('Second input argument for copyfile must be a string.')    
    end
    cmnd = ['copyfile ' source ' ' dest];
    eval(cmnd);
end