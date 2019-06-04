% [x,X_] = createAReqnMatrix(X, p, [noConstant], [future])
% If the user input a time series vector then we create the 
% x & X_ in x = X_*A + e
% Specifically, assuming t samples:
%  |X(p+1)|   |1  X(p) ...    X(1)  | | a(0) |    |e(p+1)|
%  |   :  |   |    :    ...    :    | | a(1) |    |      | 
%  |X(t-1)| = |1 X(t-2) ... X(t-p-1)| |  :   |  + |e(t-1)|
%  | X(t) |   |1 X(t-1) ... X(t-p)  | | a(p) |    | e(t) |
% (If noConstant==1 the column vector of 1's disappears).
% So this assumes that the time series is modelled by an AR model of order
% p, and we want to find the RHS vector and middle matrix to then solve this
% matrix equation by finding A, using getARcoeff.
% Inputs:
%   X   - X is a vector of time series samples
%   p   - p is the autoregression order
%   noConstant - set to 1 if you want x(t) = a1.x(t-1) + ... + ap.x(t-P), 
%         instead of x(t) = c + a1.x(t-1) + ... + ap.x(t-P)
%   future - set to 1 to make projection matrix from p current 
%            & future values
% Outputs:
%   x   - original time series without the first p entries
%   X_  - matrix form of time series with p+1 entries in each row
% See also getARcoeff, getGaussAR
function [Y,X_] = createAReqnMatrix(X, p, varargin)
    nargs = length(varargin);
    if nargs>2
        error('createAReqnMatrix:TooManyInputs', ...
            'requires at most 2 optional inputs');
    end
    optargs = {0 0}; % {noConstant}
    optargs(1:nargs) = varargin;
    [noConstant,future] = optargs{:};

    X = squeeze(X);
    if ~isvector(X) || ~isscalar(p)
        error('X must be a vector & p must be a scalar');
    end
    X=X(:); % force to column vector
    t = length(X); % length gets longest dim, which is number of samples
    if t<2*p
       error('Time series must be longer than 2*p (T=%d, p=%d)',t,p); 
    end
    if ~future
        pinv = p:-1:1; % [p p-1 p-2 ... 1]
        Y = X(p+1:end); % Y = X_*A + e
        X_ = [ones(t-p,~noConstant) zeros(t-p,p)];
        for i=1:1:p
            X_(:,~noConstant+pinv(i)) = X(i:t-pinv(i));
        end
    else
        Y = X(1:t-p+1);
        X_ = [ones(t-p+1,~noConstant) zeros(t-p+1,p)];
        for i=1:p
           X_(:,~noConstant+i) = X(i:t-p+i);
        end
    end
end
