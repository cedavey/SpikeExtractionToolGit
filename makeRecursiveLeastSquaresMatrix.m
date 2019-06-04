%   X = makeRecursiveLeastSquaresMatrix( x, y, n, na, nb, [delay] )
%
% Construct sample matrix for regression from input & output vectors for
% the system
%
%         z^-d * (bo + b1*z^-1 + b2*z^-2 + ... + b_nb*z^-nb)
% G(z) =  ---------------------------------------------------
%                1 + a1*z^-1 + ... + a_na*z^-na
%
% Inputs:
% x     - input (must be at least max(nb, na+1) long)
% y     - observations / output (must be at least max(nb, na+1) long)
%       - x & y have the most recent sample last, but when you're
%         constructing the X matrix we're doing it in the opposite order so 
%            X = [ y(n-1) ... y(n-na) x(n) ... x(n-nb+1) ]
% n     - current timestep, so can input whatever size vector you like
% na    - get na of the previous y samples (ignore current sample as we only
%         include lagged autoregression)
% nb    - get nb of the current & previous x samples (include current input)
%         (includes dc component, b0, so always at least 1)
% delay - delay (if provided) is applied to the input vector
%
% See also:    recursiveLeastSquares 
function X = makeRecursiveLeastSquaresMatrix( x, y, n, na, nb, delay )
   if nargin<6 || isempty(delay), delay=0; end
   nu = na + nb;
   
   X = zeros(nu-1,1); % -1 because we don't want DC bit yet 
   for j = 1:(nu-1)
      if j <= na % terms of y
         if (n-j)<=0
             X(j) = 0;
         else
             X(j) = y(n-j);
         end
      else       % terms of x
         if (n-delay-(j-(na+1)))<=0
             X(j) = 0;
         else
             X(j) = x(n-delay-(j-(na+1)));
         end
      end
   end
   
   X = [1; X]; % regression samples is [DC auto_regress input_regress]

end