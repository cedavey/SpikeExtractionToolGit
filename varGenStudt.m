%     v = varGenStudt( x )
% Returns the variance of samples from a generalised students-t
% distribution, such that X = mu + T*sigma. Optionally returns the
% structure containing the fit parameters for studt, from fitdist.
function [ v, t ] = varGenStudt( x, t )
   if nargin==0
      help varGenStudt;
      return;
   end
   % For generalised students-t(mu, sigma, nu) we know that 
   % var(x) = sigma^2 * nu / (nu - 2)
%    v = t.sigma^2 * t.nu / (t.nu - 2);
   if nargin<2
      t = fitdist( x, 'tlocationscale');
   end
   
   v = t.sigma ^ 2;
end