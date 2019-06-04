% zn = normStudt( x ) 
% Normalise samples from a generalised students-t distribution to be
% standard normal
function z = normStudt( x, t ) 
   if nargin==0
      help normStudt;
      return;
   end
   if nargin<2
      t  = fitdist(x,'tlocationscale'); 
   end
   zt = (x - t.mu) / t.sigma; % transform to a 1-parameter students-t distribution
   
   % Use the transformation proposed by Walcott, discussed by Prescott in 
   % unimelb -> articles -> maths -> distributions -> Prescott74_Biometrika_NormalisingStudtDist
   z  = (8*(t.nu)+1)/(8*(t.nu)+3) * ( (t.nu) * log( 1 + zt.^2/(t.nu) ) ).^0.5 .* sign(zt);
end