% rho = spikeCorrWithTemplates( matchfn, sp, templates, Ntemplates )
%
% Calculate correlation of a spike with the family spike templates.
% Templates can be a matrix if they're all forced to be the same length,
% else it's a cell array if they're allowed to differ in length
% Inputs:
%  sp         - spike to be correlated with templates
%  templates  - matrix (time x family) or cell array of vectors
%  alignment_ind_sp - index to align on for the spike 
%  alignment_ind_templates - vector of indices to align on for the
%  templates
%  Ntemplates - number of current templates
%  matchfn    - match function to calc similarity of spike with templates
function rho = spikeCorrWithTemplates( sp, templates, alignment_ind_sp,alignment_ind_templates,matchfn )
   if nargin<5 || isempty( matchfn )
      matchfn = @( sp, temp ) corr( temp, sp );
   end
   
   % templates & spikes have different lengths
   if iscell( templates )
      try
         Ntemplates = length(templates);
         rho = zeros(Ntemplates,1); %preallocate for speed
         for si = 1:Ntemplates
            [tmp,sp_] = alignDiffLengthSpikes( templates{si}, sp,[], [],alignment_ind_templates(si),alignment_ind_sp,  false );
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
       try
        Ntemplates = size(templates,2);
        rho = zeros(Ntemplates,1); %preallocate for speed
        for si = 1:Ntemplates
            rho(si,1) = matchfn( sp, templates(:,si));
        end
       catch ME
         str = getCatchMEstring( ME );
         cprintf( 'keyword*', str );
         runtimeErrorHandler( ME );
         return;
      end
   end
end











