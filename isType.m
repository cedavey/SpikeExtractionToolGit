% isType(var, [type])
% This function determines the type of the input variable. Can be used in
% one of 2 different modes. If the user specifies a type in the input then	
% the function returns true if the variable matches the specified type, or
% false otherwise. If no type is specified then the function returns the
% variable type in a string.
% The different types checked for:
%   empty     - checks if the input is empty
%   char      - a single character
%   scalar    - a single number (integer, floating point or logical)
%   string    - a vector of characters
%   vector    - a vector of numbers (integer, floating point Or logical)
%   matrix    - a matrix of numbers (integer, floating point or logical)
%   structure - a structure containing any type of elements
function istype = isType(var, type)
    istype=[];
    if nargin<2
        if isempty(var)
    %       istype = 'empty'
        elseif isscalar(var) && ischar(var)
            istype = 'char';
        elseif isvector(var) && ischar(var)
            istype = 'string';
        elseif isfloat(var)
            istype = 'float';
        elseif islogical(var)
            istype = 'boolean';
        elseif isinteger(var)
            istype = 'int';
        elseif isscalar(var) && (isfloat(var) || isinteger(var) || islogical(var))
            istype = 'scalar';
        elseif isvector(var) && (isfloat(var) || isinteger(var) || islogical(var))
            istype = 'vector';
        elseif (isfloat(var) || isinteger(var) || islogical(var))
            % if here then float/int matrix, since not scalar or vector
            istype = 'matrix'; 
        elseif isstruct(var)
            istype = 'struct';
        elseif iscell(var)
            istype = 'cell';
        elseif isstruct(var)
            istype = 'struct';
        else
            istype = 'other';
        end
        return;
    else % user provided type to check for
        istype = 0;
        switch lower(type)
            case {'empty'}
                if isempty(var)
                    istype = 1;
                end
            case {'char','character'}
                if isscalar(var) && ischar(var)
                    istype = 1;
                end
            case {'float','double'}
                if isfloat(var)
                    istype = 1;
                end
            case {'numeric'}
                if isnumeric(var)
                    istype = 1;
                end
            case {'integer','int'}
                if isint(var)
                    istype = 1;
                end
            case {'string','str','strings'}
                if isvector(var) && ischar(var)
                    istype = 1;
                end
            case {'scalar'}
                if isscalar(var) && (isfloat(var) || isinteger(var))
                    istype = 1;
                end
            case {'matrix','mat'}
                if (~isempty(var) && ~(any(size(var)==1))) && (isfloat(var) || isinteger(var))
                    istype = 1;
                end
            case {'vector','vec'}
                if isvector(var) && (isfloat(var) || isinteger(var))
                    istype = 1;
                end
            case {'structure','struct'}
                istype=isstruct(var);
            case {'boolean','bool','logical','logic','log'}
                istype=islogical(var);
            otherwise
                istype = 0;
        end
    end % end if type is provided by user or not
end