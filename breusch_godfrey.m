% [H, pval] = breusch_godfrey(X, [p], [alpha])
% The Breusch Godfrey test just applies an AR model to the data, then
% calculates the multiple correlation coefficient, and tests it with a chi
% squared test for significance (NR^2 ~ X^2(p)).
% H0: data is random (not serially correlated)
% H1: data contains serial correlation (is not random)
% Inputs:
%   X     - time series vector or matrix of dimensions x samples
%           if X is a matrix calc Q stat for each row
%   p     - number of lags to include in test scalar (default=3)
%   alpha - confidence interval (e.g. alpha = 0.05)
% Outputs:
%   H     - 0 for null hypothesis accepted, 1 for alternative hyp. 
%   pval  - probabilit(ies) of the chi-squared distribution

% I'm not sure why you'd test R^2 as a chi-square statistic? Perhaps I'll 
% try a few different distributions - e.g. those in my Granger Danger 
% article, & see which gives the best ROC curve. I guess the R^2 stat in my
% article, tested as an F stat, has an R^2 var in numerator & denominator -
% since this has only a single R^2, it's X^2 instead of F?
function [H,pval] = breusch_godfrey(X, varargin)
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
    
    % now we have default/user set options, all initialised, we can start
    T       = size(X,2);        % length of time series
    pval    = zeros(size(X,1),1);
    H       = zeros(size(X,1),1);
    
    
    % calculate stat for each row
    for i=1:size(X,1)
        % get lagged time series of X & Y
        [x,X_]  = createAReqnMatrix(X(i,:), p, 1);
        [R, R2] = multCorr(x, X_);      % R_x|x_{t-p}
%         pval(i) = fcdf((T-p)/p*R2/(1-R2),p,T-p);  % NR^2 ~ X^2 (N=dof)
        pval(i) = chi2cdf((T-p)*R2,p);  % NR^2 ~ X^2 (N=dof)
        
        H(i) = ((1-pval(i))<alpha);
        H(i) = ~H(i);
    end    
return % end of function