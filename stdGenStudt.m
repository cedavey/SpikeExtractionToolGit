%     [v, t] = stdGenStudt( x )
% Returns the sigma of samples from a generalised students-t
% distribution, such that X = mu + T*sigma. Optionally returns the
% structure containing the fit parameters for studt, from fitdist.
function [sig, t] = stdGenStudt( x )
   if nargin==0
      help varGenStudt;
      return;
   end
   
   x(isnan(x)) = []; % remove NaNs so it doesn't derail the fit!
   
   t = fitdist( x, 'tlocationscale');
   sig = t.sigma;

%    sig = t.nu;
end
