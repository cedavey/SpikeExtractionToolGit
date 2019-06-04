% [ok,res] = checkStringInput(str, type, [min], [max], [def])
% def is default value, wh must be a string
function [ok,res] = checkStringInput(str,type,varargin)
    nargs = length(varargin);
    optargs = {[],[],[]};
    optargs(1:nargs) = varargin;
    [min,max,def] = optargs{:};
    type = lower(type);
    
    ok = 0; res = [];
    if isempty(str)
        if isempty(def)
            return
        else
            str = def;
        end
    end
    if strcmpi(type,'char')
        ok = 1; res = str;
        return
    end
    num = str2num(str);
    if isempty(num) % empty if string is not numeric
        return
    end
    if ~isempty(min) && num<min
        return
    end
    if ~isempty(max) && num>max
        return
    end
    if strcmpi(type,'float') || strcmpi(type,'double')
        ok = 1; res = str2num(str);
        return
    end
    % not string, float or double, so must be integer or boolean, which
    % means the number must round to itself
    if round(num)~=num
        return
    end
    ok = 1; res = str2double(str);
return
