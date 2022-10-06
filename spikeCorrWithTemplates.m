% rho = spikeCorrWithTemplates( matchfn, sp, templates, Ntemplates )
%
% Calculate correlation of a spike with the family spike templates.
% Templates can be a matrix if they're all forced to be the same length,
% else it's a cell array if they're allowed to differ in length
% Inputs:
%  sp         - spike to be correlated with templates
%  templates  - matrix (time x family) or cell array of vectors
%  Ntemplates - number of current templates
%  matchfn    - match function to calc similarity of spike with templates
function rho = spikeCorrWithTemplates( sp, templates, matchfn )
   if nargin<3 || isempty( matchfn )
      matchfn = @( sp, temp ) corr( temp, sp );
   end
   
   % templates & spikes have different lengths
   if iscell( templates )
      try
         Ntemplates = length(templates);
         for si = 1:Ntemplates
            [tmp,sp_] = alignDiffLengthSpikes( templates{si}, sp, [], [], false );
            rho(si,1) = matchfn(tmp, sp_);
         end
         
      catch ME
         str = getCatchMEstring( ME );
         cprintf( 'keyword*', str );
         runtimeErrorHandler( ME );
         return;
      end
      
   % templates & spikes are all the same length
   else
      rho = matchfn( sp, templates );

   end
end











