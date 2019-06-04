% str = title2Str(str,noDecimal,noPunct,replace)
% Convert a title string or cell of strings to a string with no white spaces
% Set noDecimal to 1 when you want to convert decimal points to _ also
% Set noPunct to 1 if everything other than alphabet chars (i.e. all
% punctuation chars) should be set to _, or 'replace', if provided
% Inputs:
%   str       - string to modify
%   noDecimal - convert all decimal points to replacement char
%   noPunct   - convert all punctuation to replacement char
%   replace   - character to convert all decimal, punctuation etc to
function str = title2Str(str, varargin)
    nargs = length(varargin);
    if nargs>3
        error('title2Str:tooManyInputs','max of 3 optional args');
    end
    optargs = {1 0 '_'};
    optargs(1:nargs) = varargin;
    [noDecimal, noPunct, replace] = optargs{:};
    if isempty(replace)
        replace = ' ';
        remove = 1;
    else 
        remove = 0;
    end
    if iscell(str)
        n = prod(size(str));
        for i=1:n
            str{i}(isspace(str{i}))=replace;
            if noDecimal==1
                str{i}(isspace(str{i}))=replace;
            end
            if noPunct==1
                str{i}(~isstrprop(str,'alphanum'))=replace;
            end
        end
    else
        str(isspace(str))=replace;
        if noDecimal==1
            str(isspace(str))=replace;
        end
        if noPunct==1
            str(~isstrprop(str,'alphanum'))=replace;
        end
    end
    if remove
        str = str(str~=' ');
    end
return