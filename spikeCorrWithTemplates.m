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
   if nargin<4 || isempty( matchfn )
      matchfn = @( sp, temp ) corr( temp, sp );
   end
   
   % templates & spikes have different lengths
   if iscell( templates )
      try
         % if spikes & templates can have diff lengths, use the smaller length
         % to calculate the match, but centre both around the peak 
         Ntemplates = length( templates );
         nSp     = length( sp );
         nT      = cellfun( @length, templates ); % length of each template
         peakind = cellfun( @getMaxInd, templates ); 
         peaksp  = getMaxInd( sp );
         tbefore = min( [ peakind(:) ones(Ntemplates,1)*peaksp ], [], 2 );
         tafter  = min( [ nT(:) - peakind(:), ones(Ntemplates,1)*(nSp - peaksp) ], [], 2 );

         % if correlating on less than half the template then assume something is
         % wacky with the spike
         if max( tbefore + tafter ) < min( cellfun( @length, templates ) ) / 2 ...
         || all( tafter == 0 ) || all( tbefore == 0 )
            rho = [];
            return;
         end

         % gotta calc corr or cov separately for each template because
         % they're different lengths
         rho = zeros( Ntemplates, 1 );
         for si = 1:Ntemplates
            rho(si) = matchfn( sp( peaksp-tbefore(si)+1 : peaksp+tafter(si) ), ...
                               templates{si}( peakind(si)-tbefore(si)+1 : peakind(si)+tafter(si) ) );
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











