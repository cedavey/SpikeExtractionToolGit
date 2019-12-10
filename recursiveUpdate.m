% [mu_curr, var_curr] = recursiveUpdate( mu_prev, var_prev, lambda, new_obs )
% 
% OR
% 
% mu_curr = recursiveUpdate( mu_prev, [], lambda, new_obs )
%
% Recursive online calculation of mean and variance, with optional
% forgetting factor, lambda. Lambda is typically between 0.9 & 0.99, where
% the larger lambda is the more slowly the mean and variance are allowed to
% change (the smoother the estimates change over time)
function [mu_curr, var_curr] = recursiveUpdate( mu_prev, var_prev, new_obs, lambda )
   if nargin<5 || isempty( lambda ), lamba = 1; end
   
   % mu_curr  = mu_prev * lambda  +  (mu_prev + (new_obs - mu_prev)) * (1-lambda);
   
   % https://math.stackexchange.com/questions/374881/recursive-formula-for-variance
   mu_curr  = mu_prev  +  (new_obs - mu_prev) * (1-lambda);
   
   if ~isempty( var_prev )
      % these 2 var estimates give the same result
      % var_curr = var_prev + ( mu_prev^2 - mu_curr^2 ) + ( new_obs^2 - var_prev - mu_prev^2 ) * (1-lambda);
      var_curr = var_prev * lambda + ( new_obs - mu_prev ) * ( new_obs - mu_curr ) * (1-lambda);
   end
   
end

% testing
% lambda=0.9; mux=x(1); sigx=1; 
% for i=2:length(x)
%    [mux(i), sigx(i)] = recursiveUpdate( mux(i-1), sigx(i-1), x(i), lambda ); 
% end; 
% figure; hold on; errorbar( mux, sqrt(sigx) ); plot(x,'linewidth',1); 
% [mean(x) mean(mux) var(x) mean(sigx)]