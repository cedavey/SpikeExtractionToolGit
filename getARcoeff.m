% [A, e, e2] = getARcoeff(X, p, [x], [doplot], [warn])
% Returns the least squares estimate of an AR time series (xt), of order
% p. Finds the A matrix in:
%    x(t) = a(1)*x(t-1) + ... + a(p)*x(t-p) + e, where e ~ N(0,sigma)
% The user can input the data in 2 ways - either they can provide a 
% time series vector:     [X(1), X(2), ... X(t-1), X(t)] 
% or they can provide a matrix of equations. So, assuming t samples:
%  |x(p+1)|   |1  x(p) ...    x(1)  | |a(p)|    |e(p+1)|
%  |   :  |   |    :    ...    :    | | :  |    |      | 
%  |x(t-1)| = |1 x(t-2) ... x(t-p-1)| |a(2)|  + |e(t-1)|
%  | x(t) |   |1 x(t-1) ... x(t-p)  | |a(1)|    | e(t) |
%                                     |a(0)|
% which can be written more generally as x = XA + e the user can provide
% Y and X. Both input methods require the user to input the order, p.
% (Perhaps a little confusing having 2 input formats but it's designed to
% avoid creating the matrix every time if using the function iteratively to 
% determine p).
% You can choose to throw a warning if X is singular. If you want to throw
% a warning but not provide the Y matrix, make Y an empty matrix. 
% Returns all zeros if X is singular and A can't be found.
% Input:
%   X     - either a vector of samples or a matrix of eqns (the X in x = XA + e)
%   p     - order of the AR. p coefficients will be returned.
%   x     - the x in x = XA + e, only required if inputting a matrix of eqns
%           instead of a vector of samples.
%   warn  - if 1 throw a warning if X matrix is (almost) singular (default=0)
% Output:
%   A     - the p+1 coefficients estimated using least squares (p0
%           represents a shift from the mean
%   e     - (vector) the residual error for each equation given A
%	e2    - average error squared (i.e. sum of errors squared/num samples
function [A, e, e2] = getARcoeff(X, p, varargin)
    if nargin==0
        help getARcoeff;
        return;
    end
    optargs = {[],0, 0}; % {Y, doplot, warn}
    na = length(varargin);
    if na>3
        error('getARcoeff:TooManyInputs','requires at most 2 optional args');
    end
    optargs(1:na)   = varargin;
    [Y,doplot,warn] = optargs{:};
    
    if ~isscalar(p) || p<0
        error('order of AR function (p) must be a non-negative scalar value');
        return
    end
    % if AR(0) then error is original time series
    if p==0
        e = X;
        A = 0;
        e2 = e'*e/(length(e)-1);
        return
    end
    X = squeeze(X);   % if got 1 pixel out of a volume if may have dim 1x1x1x200
    if isempty(Y) && ~isvector(X) 
        error('if inputting X & p only then X must be a vector of samples');
        return
    end
    if isinteger(X)   % X gotta be type double to perform matrix ops
        X = double(uint16(X));
    end
    if isempty(Y) % input is a time series vector, X
        t      = length(X); % length gets longest dim of X - i.e. num samples
        [Y, X] = createAReqnMatrix(X, p);
    elseif isscalar(Y) % I keep putting doplot where Y is, so just switch if I do
        doplot = Y;
        t      = length(X); % length gets longest dim of X - i.e. num samples
        [Y, X] = createAReqnMatrix(X, p);
    else
        if ~( isDim(X,2) || ~isvector(Y) )
            error('error providing matrix input - see help for details')
            return
        end
        t = p + size(X,1);
        if size(X,1) ~= size(Y,1)
            error('X & Y must be column vectors of equal length - see help for details')
            return
        end
        if size(X,2)~=p+1
            error('X must have p+1 columns - see help for details')
            return
        end
    end
    if ~warn
        warning off MATLAB:singularMatrix;
        warning off MATLAB:nearlySingularMatrix;
        warning off MATLAB:divideByZero;
    end
    
    try % x = XA + e --> A = X\x
        % note: this is equivalent to A = inv(X'*X)*X'*Y
        A = X\Y;  % A = (X'X)^-1*X'Y (from min. squared error)
        if isinf(A) | isnan(A) % throw error manually cuz' removed warnings
            A = 0; e = 0; e2 = 0; % sometimes the catch isn't thrown from 'error'??
            error('');
        else
            e  = Y - X*A; 
            e2 = e'*e/(length(e)-1);
        end
    catch
        A = 0;
        e = 0;
        e2 = 0;
        if warn
           warning('Matrix is singular - cannot estimate coefficients');
        end
    end
    
    if doplot
        figure; 
        plot(e,'r--','LineWidth',1);
        hold on; 
        plot([Y X*A]); 
        legend({'error','x','x predict'});
    end
    
    if ~warn
        warning on MATLAB:singularMatrix;
        warning on MATLAB:nearlySingularMatrix;
        warning on MATLAB:divideByZero;
    end
end
