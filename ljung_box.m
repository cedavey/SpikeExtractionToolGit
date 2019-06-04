% [H, pval, qstat] = ljung_box(X, [p], [alpha])
% Computes the Ljung-Box Q test for lags up to & including AR(p). 
% H0: data are random (uncorrelated)
% H1: data contains serial correlation (not random)
% Inputs:
%   X     - time series vector (usually regression residuals)
%           if X is a matrix calc Q stat for each row
%   p     - number of lags to include in test (scalar or vector) (default=10)
%   alpha - confidence interval (e.g. alpha = 0.05)
% Outputs:
%   H - set to 1 if satisfies null hypothesis (chi-squared dist) 
%   qstat - Q(p) test statistic, distributed Chi-squared(p) 
%           under null hypothesis. Gives the test statistic value for each
%           lag up to the max lag (in increasing order)
%   pval  - probabilit(ies) of the qstat(s)
function [H, pval, qstat] = ljung_box(X, varargin)
    % only want 2 optional inputs at most
    numvarargs = length(varargin);
    if numvarargs > 2
        error('ljung_box:TooManyInputs', ...
            'requires at most 2 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {3 0.05}; % optargs = {p alpha}
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite with the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    % [optargs{1:numvarargs}] = varargin{:}; % or this format
    % Place optional args in memorable variable names
    [p, alpha] = optargs{:};
    
    % if X is vector force to row vector - else many dims with sample size 1
    if isvector(X)
        X = X(:)';
    end
    if isempty(p)
        p=3;
    end
    
    % now we have default/user set options, all initialised, we can start
    p     = sort(p(:));       % force p to be a column vector
    T     = size(X,2);        % length of time series
    n     = size(X,1);
    qstat = zeros(n,length(p));
    pval  = zeros(n,length(p));
    H     = zeros(n,1);
    
    
    % calculate stat for each row
    for i=1:n
        rho   = acf(X(i,:), p(end));   % calc ACF
        rho   = rho(p(1)+1:end); % ignore ACF(0) cuz' nec equals 1
        qstat(i,:) = T.*(T+2).*cumsum(rho.^2/(T-p));  % calc Q (ljung-box stat)
        pval(i,:)  = chi2cdf(qstat(i,:)',p);          % calc prob of Q

        H(i) = ((1-pval(i,end))<alpha);
        H(i) = ~H(i);
    end    
return % end of functionx

