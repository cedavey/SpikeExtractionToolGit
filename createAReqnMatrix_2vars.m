% [z, ZT] = createAReqnMatrix_2vars(X, Y, p)
% Creates matrices out of time series in order to solve a bivariate AR
% problem. X & Y are 2 input time series (vectors) that were created from
% (or are being modeled by):
%  x(k) = a(0) + a(1)*x(k-1) + ... + a(p)*x(k-p) + b(1)*y(k-1) + ... + b(p)*y(k-p) + e(k)
%  y(k) = c(0) + c(1)*y(k-1) + ... + c(p)*y(k-p) + d(1)*x(k-1) + ... + d(p)*x(k-p) + w(k)
% If X & Y are different lengths the longest is shortened by ignoring the
% end values, since for each equation we need samples from both time
% series. So both time series then have T samples.
% To find the coefficients {a,b,c,d} by least squares, the following matrix
% equation is constructed (where a is vector form of {a(1),...,a(p+1)}), etc):
% |X(p+1)|   |  1    X(p) ... X(1) 0(p+1) 0(p) ... 0(1)   Y(p) ... Y(1)  0(p) ... 0(1) | |a|   |e(p+1)|
% |Y(p+1)|   |0(p+1) 0(p) ... 0(1)   1    Y(p) ... Y(1)   0(p) ... 0(1)  X(p) ... X(1) | |c|   |w(p+1)|
% |   :  | = |  :     :   ...   :    :     :   ...  :      :   ...  :     :   ...   :  | |b| + |   :  |
% | X(T) |   |  1   X(T-1)...X(T-p) 0(p+1) 0(p)... 0(1)  Y(T-1)...Y(T-p) 0(p) ... 0(1) | |d|   | e(T) |
% | Y(T) |   |0(p+1) 0(p) ... 0(1)   1   Y(T-1)...Y(T-p)  0(p) ... 0(1) X(T-1)...X(T-p)|       | w(T) |
%
% If the 2 AR equations have different orders then we get the following
% (note that x & y within each eqn are still assumed to have the same
% order). W.l.o.g assume order of y(k) eqn is larger than x(k) eqn for the e.g.
%|X(px+1)|   |  1    X(px) ... X(1) 0(py+1) 0(py)  ... 0(1)  Y(px)... Y(1)  0(px) ...  0(1) | |a(0) |   |e(px+1)|
%|   0   |   |  0      0   ...  0      0      0    ...   0    0   ...   0     0   ...   0   | |  :  |   |   0   |
%|   :   |   |  :      :   ...  :      :      :    ...   :    :   ...   :     :   ...   :   | |a(px)|   |   :   |
%|X(py+1)| = |  1    X(px) ... X(1) 0(py+1) 0(py)  ... 0(1)  Y(px)... Y(1)  0(py) ...  0(1) | |c(0)|    |e(py+1)|
%|Y(py+1)|   |0(px+1) 0(px)... 0(1)    1    Y(py)  ... Y(1)  0(px)... 0(1)  X(py) ...  X(1) | |  :  |   |w(py+1)|
%|   :   |   |  :      :   ...  :      :      :    ...   :    :   ...   :     :   ...   :   | |c(py)| + |   :   |
%| X(T)  |   |  1    X(T-1)...X(T-px) 0(py+1) 0(py)... 0(1) Y(T-1)...Y(T-px) 0(py)...  0(1) | |b(0) |   | e(T)  |
%| Y(T)  |   |0(px+1) 0(px)... 0(1)    1    Y(T-1)...Y(T-py) 0(px)... 0(1)  X(T-1)...X(T-py)| |  :  |   | w(T)  |
%                                                                                             |b(px)|
%                                                                                             |d(0) |
%                                                                                             |d(py)|
% In matrix notation: z = Z'A + E = ZT*A + E (z, A, E are vectors, Z a matrix)
% Inputs:
%   X   - (vector) periodic samples of the AR process
%   Y   - (vector) periodic samples of the AR process
%   p   - (scalar) AR order for which to find matching coefficients
% Outputs:
%   z   - LHS time series vector that alternates btwn x & y samples
%   ZT  - the big RHS matrix of samples from X & Y

% FIXME: allow mean shift / assume it's been removed / remove it ???
% FIXME: allow x(k) & y(k) eqns to have different orders (although x & y
% within each eqn gotta have same order - coeffs just set to 0 if unused)??
function [z,ZT] = createAReqnMatrix_2vars(X, Y, p)
    X = squeeze(X); Y = squeeze(Y);
    if ~isvector(X) | ~isvector(Y) | ~isscalar(p)
        error('X & Y must be vectors, p must be scalar');
    end
    if size(X,1) < size(X,2), X = X'; end % convert to column vector
    if size(Y,1) < size(Y,2), X = X'; end % convert to column vector
    if length(X)<=length(Y), T = length(X); Y = Y(1:T); end
    if length(X)> length(Y), T = length(Y); X = X(1:T); end
    
%    p = max(px, py); 
    if T < (4*p+2)
        error('length of time series must be greater than 4*p+2 for non-singular matrix');
    end
    
    % create the LHS vector, z, from z = Z'A + e
    z = zeros(2*(T-p),1); % x elmt & y elmt in each of T-p equations
    z(1:2:end-1) = X(p+1,:); % Z alternates btwn x & y elmts
    z(2:2:end)   = Y(p+1,:);
    
    % create the RHS matrix, Z, from z = Z'A + e
    ZT    = zeros(2*(T-p), 4*p+2);
    vec0  = zeros(p,1); % for the runs of 0's in each Z row
    for i = 1:(T-p) % each time series has T-p equations
        t = i + p; % sample = iteration + p (cuz eqns start from order+1)
        z(2*i-1)    = X(t);
        z(2*i  )    = Y(t);
% size(ZT(2*i-1,:))
% size([1; X(t-1:-1:t-p)])
% size(vec0)
% size([1; Y(t-1:-1:t-p)])
% size(vec0)
        ZT(2*i-1,:) = [1; X(t-1:-1:t-p);    0; vec0;       Y(t-1:-1:t-p);    vec0      ];
        ZT(2*i  ,:) = [   0; vec0;       1; Y(t-1:-1:t-p);     vec0;      X(t-1:-1:t-p)];
    end
end
