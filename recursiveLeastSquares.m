%  [a, b, P, y_est, g] = recursiveLeastSquares(x, y, na, nb, P_prev, a_prev, b_prev, lambda, n, delay)
% 
% Recursively calculate the coefficients for a known system, defined by the 
% following transfer function:
%
%         z^-d * (bo + b1*z^-1 + b2*z^-2 + ... + b_nb*z^-nb)
% G(z) =  ---------------------------------------------------
%                1 - a1*z^-1 - ... - a_na*z^-na
%
% Inputs:
%     x   - input data vector
%     y   - output vector (observations)
%     na  - number of poles in the transfer function (all lagged, so that
%           obs vector must be at least na+1 samples long)
%     nb  - number of zeros in the transfer function
%     P   - weighted inverse covariance matrix
%     a   - autoregression coeffs that minimise the least squares cost fn
%     b   - DC + input coeffs that minimise the least squares cost function
% lambda  - forgetting factor; smaller lambda --> smaller the contribution
%           of previous samples to the covariance matrix
%     n   - current timestep index, so we only use 'delay + nb' previous 
%           input values and 'na' output values
%     d   - is the delay on the input 
% Outputs:
%     a   - estimate of autoregressive coefficients (i.e. coeffs of lagged
%           output)
%     b   - estimate of regression coefficient (i.e. coeffs of input)
%     P   - updated inverse covariance matrix
%     y_est - estimated value of observation y using coeffs a & b
%     w   - updated coefficient vector, from which a & b are drawn
%           such that w = [1 a(1) a(2) ... a(na) b(1) b(2) ... b(nb)]
%           which is regressed on [dc y(n) y(n-1) ... y(n-nb) x(n) x(n-1) ... x(n-na)]
%     g   - the gain vector, so that if alpha is the error between the
%           output & our estimate of the output, then the new coefficient
%           estimate is w <- w + alpha*g
%     X   - updated matrix created from input & output vectors, for regression
% a and b are vectors of the system estimated parameters.
%
function [a, b, P, y_est, g] = recursiveLeastSquares(x, y, na, nb, P_prev, a_prev, b_prev, varargin)
   optargs = {1, [], 0, 1}; % lambda, n, delay
   nargin  = length(varargin);
   optargs(1:nargin) = varargin(:);
  [lambda, n, delay] = optargs{:};
   nu = na + nb; % number of unknowns that we need to estimate (dc is in nb)

   if isempty(lambda), lambda = 1; end
   if isempty(n),   n = max(length(x),length(y)); end
   if isempty(delay),   delay = 0; end

   if length(y) < max(nb, na+1)
      str = sprintf('\tOutput needs at least nb, or na+1 (%d) samples\n', na);
      cprintf('Error', str);
      a=a_prev; b=b_prev; P=P_prev; y_est=[]; g=[];
      return;
   end
   if n < (delay+nb)
      str = sprintf('Output needs at least delay + nb (%d) samples\n', delay+nb);
      cprintf('Keywords*', str);
      a=a_prev; b=b_prev; P=P_prev; y_est=[]; g=[];
      return;
   end
   if length(P_prev) ~= nu   ||   size(P_prev,1) ~= size(P_prev,2)
      str = sprintf( 'Inverse cov matrix P must be square with dimension (na + nb)\n' );
      cprintf('Keywords*', str);
      a=a_prev; b=b_prev; P=P_prev; y_est=[]; g=[];
      return;
   end

   X = makeRecursiveLeastSquaresMatrix( x, y, n, na, nb, delay );

   % create prev weights to match regression matrix: [ DC auto_regress input_regress ] 
   w_prev = [b_prev(1); a_prev(:); toVec(b_prev(2:end))]; 
   
   %% Estimation
   % https://en.wikipedia.org/wiki/Recursive_least_squares_filter#LRLS_Algorithm_Summary
   alpha = y(end) - X'*w_prev;
   g = P_prev * X * inv( lambda + X'*P_prev*X );
   P = lambda^(-1) * P_prev   -  (g) * X' * lambda^(-1) * P_prev;
   w = w_prev + alpha * g;
   % System Parameters
   a =  w(2:na+1); % w(1:na);
   b = [w(1); w(na+2:end)]; % w(na+1:end)
   
   y_est = X'*w;
end



