% [H p] = iswhite(data,[test],[alpha],[p])
%
% Test each time series (or image, if spatial) for whiteness.
% Inputs
%   data - vox x samples
%   test - white test to use
%          'ljung-box' (lb), 'breusch-godfrey' (bg)
%   p    - number of lagged samples to use in AR fit
% Outputs
%   H    - whiteness hypothesis supported (1) or rejected (0)
%
% See also: isnormal, isstationary
function [H, pval] = iswhite(data, varargin)
    nargs = length(varargin);
    if nargs>3
        error('isnormal:inputError','max 2 optional args');
    end
    optargs = {'lb', 0.05, 3}; % test=anderson-darling, order=1, alpha=0.05
    optargs(1:nargs) = varargin;
    [test, alpha, p] = optargs{:};
    if isempty(test)
        test = 'bg';
    end

    if isvector(data)
        data = data(:)';
    end
    if isempty(p)
        p=3;
    end
    
    n = size(data,1);
    H = zeros(n,1);
    pval = zeros(n,1);

    switch lower(test)
        case {'ljung-box','lb','l'} % H=1: data random, H=0: not random
            [H, pval] = ljung_box( data, p, alpha );
            
        case {'breusch-godfrey','bg','b'} % H=1: data random, H1: not random
            [H, pval] = breusch_godfrey( data, p, alpha );
    end
return