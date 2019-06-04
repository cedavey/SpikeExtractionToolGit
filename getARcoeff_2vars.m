% [Ax, Ay, ex2, ey2, ex, ey] = getARcoeff_2vars(X, Y, p, [doplot], [warn])
% Returns the least squares estimate of an AR bivariate time series of order
% p. Finds the coefficient vectors a, b, c, d from:
%  x(k) = a(0) + a(1)*x(k-1) + ... + a(p)*x(k-p) + b(1)*y(k-1) + ... + b(p)*y(k-p) + e(k)
%  y(k) = c(0) + c(1)*y(k-1) + ... + c(p)*y(k-p) + d(1)*x(k-1) + ... + d(p)*x(k-p) + w(k)
% If X & Y are different lengths the longest is shortened by ignoring the
% end values, since for each equation we need samples from both time
% series. So both time series then have T samples.
% To find the coefficients {a,b,c,d} by least squares, the following matrix
% equation is constructed (where A is vector form of {a(1),...,a(p+1)}), etc):
% |X(p+1)|   |  1    X(p) ... X(1) 0(p+1) 0(p) ... 0(1)   Y(p) ... Y(1)  0(p) ... 0(1) | |a|   |e(p+1)|
% |Y(p+1)|   |0(p+1) 0(p) ... 0(1)   1    Y(p) ... Y(1)   0(p) ... 0(1)  X(p) ... X(1) | |c|   |w(p+1)|
% |   :  | = |  :     :   ...   :    :     :   ...  :      :   ...  :     :   ...   :  | |b| + |   :  |
% | X(T) |   |  1   X(T-1)...X(T-p) 0(p+1) 0(p)... 0(1)  Y(T-1)...Y(T-p) 0(p) ... 0(1) | |d|   | e(T) |
% | Y(T) |   |0(p+1) 0(p) ... 0(1)   1   Y(T-1)...Y(T-p)  0(p) ... 0(1) X(T-1)...X(T-p)|       | w(T) |
% 
% Letting A = {a,c,b,d}, z = LHS vector, Z' = RHS matrix, & n = noise
% vector, the least squares solution of z = Z'A + n is:
%                           A = inv(Z*Z')*Z*z
% Inputs:
%   X     - time series vector of periodic samples
%   Y     - time series vector of periodic samples
%   p     - proposed AR order for which to find coefficients
%  doplot - set to true to plot original data & estimate 
%  warn   - display a warning if the matrix is non-singular (can't find coeffs)
% Outputs:
%   Ax    - vector of coefficients for the x(k) equation
%   Ay    - vector of coefficients for the y(k) equation
%  ex2    - avg sum of squared error for the x(k) equation
%  ey2    - avg sum of squared error for the y(k) equation
%   ex    - vector of errors associated with the x(k) equation
%   ey    - vector of errors associated with the y(k) equation

% FIXME: createAReqnMatrix_2vars & getARcoeff_2vars currently assume that
% the order of the x(k) & y(k) equations are both p. If assumption is
% incorrect will the extra coefficients automatically set themselves to 0,
% or will the answer just be flat out wrong??
function [Ax, Ay, Xest, Yest, ex, ey] = getARcoeff_2vars(X, Y, p, doplot, warn)
    if nargin==0
        help getARcoeff_2vars;
        return;
    end
%    noY = 0; % records if Y was provided by user, or if we have to build it
    if nargin<5 || isempty(warn),     warn = 0;  end
    if nargin<4 || isempty(doplot), doplot = 0;  end
    if nargin<3
        error('usage: X, Y & p must be provided - see help for details');
    end
    if ~isscalar(p) || p<0
        error('order of AR function (p) must be a non-negative scalar value');
    end
    if length(X)<=length(Y), T = length(X); Y = Y(1:T); end
    if length(X)> length(Y), T = length(Y); X = X(1:T); end
    if T < (4*p+2)
        error('length of time series must be greater than 4*p+2 for non-singular matrix');
    end
    X = squeeze(X); % if got 1 pixel out of a 4D volume x time matrix
    if size(X,1) < size(X,2), X = X'; end % convert to column vector
    if size(Y,1) < size(Y,2), X = X'; end % convert to column vector
    
    % X & Y both gotta be type double to perform matrix ops
    if isinteger(X), X = im2double(uint16(X));  end
    if isinteger(Y), Y = im2double(uint16(Y));  end

    [z,ZT] = createAReqnMatrix_2vars(X, Y, p);
    Z = ZT';
    
    warning off MATLAB:singularMatrix;
    warning off MATLAB:nearlySingularMatrix;
    warning off MATLAB:divideByZero;

    try
        A = inv(Z*ZT)*Z*z; % equivalent to Z'\z
        if isinf(A) | isnan(A) | isempty(A) % throw error to enter catch
            % don't actually need this line cuz' set in catch but just in case
            Ax = 0; Ay = 0; ex = NaN; ey = NaN; ex2 = NaN; ey2 = NaN;
            error('force catch'); % should enter catch with this but doesn't always??
        else
            e = z - ZT*A;
            %  A = {a0,a,c0,c,b,d} (a0 & c0 are mean shift components)
            % Ax = {a0,a,b}, Ay = {c0,c,d}
            % separate out the 2 different AR equations
            Ax = [A(1:p+1);       A(2*(p+1)+1:3*p+2)];
            Ay = [A(p+2:2*(p+1)); A(3*p+3:end)      ];
            ex = e(1:2:end-1);
            ey = e(2:2:end);
            ex2 = ex'*ex/(length(ex)-1); % -1 so it's an unbiased estimator of variance
            ey2 = ey'*ey/(length(ey)-1); % (assuming error is 0-mean)

            Zest = ZT*A;
            Xest = Zest(1:2:end); 
            Yest = Zest(2:2:end); 

            if doplot
              figure; 
                subplot(211), 
                  plot( X(:) ); hold on; plot( Xest ); 
                   legend( 'X', 'X_{est}' ); % title( 'X' );
                subplot(212), 
                  plot( Y(:) ); hold on; plot( Yest ); 
                  legend( 'Y', 'Y_{est}' ); % title( 'Y' ); 
            end
        end
    catch
        Ax = 0; Ay = 0; ex = NaN; ey = NaN; ex2 = NaN; ey2 = NaN;
        if warn
           warning('Matrix is singular - cannot estimate coefficients');
        end
    end
    warning on MATLAB:singularMatrix;
    warning on MATLAB:nearlySingularMatrix;
    warning on MATLAB:divideByZero;
end
