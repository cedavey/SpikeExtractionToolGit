% sp = alignSpikeToTemplate( temp, sp )
%
% Align a spike to a template so that the peaks are aligned, and the spike
% is the same length as the template (with NaNs added if necessary)
%
% Inputs:
%    temp  - vector template
%      sp  - vector spike
function sp = alignSpikeToTemplate( temp, sp )

   peaktemp  = getMaxInd( temp ); % find index of template peak 
   peakspike = getMaxInd( sp   ); % find index of spike peak 

   nSp      = length( sp   );    % length of spike
   nT       = length( temp );    % length of template
   tbefore  = max( peaktemp, peakspike );
   tafter   = max( nT - peaktemp, nSp - peakspike );

   % make spike match template by aligning its peak with the templates
   if peaktemp < peakspike
      % remove initial 
      sp(1:(peakspike - peakind)) = [];
   end
   
   % append NaNs to template spikes or new spike to make them the same size
   prependfn= @( sp, N ) [ NaN( N, size(sp,2)); sp ];
   appendfn = @( sp, N ) [ sp; NaN( N, size(sp,2)) ];
   % spike peak occurs before template peak
   if peakspike < tbefore 
      % if spike has fewer samples before peak, prepend NaNs 
      sp = prependfn( sp, tbefore - peakspike );
   elseif peaktemp < tbefore
   end
   if nSp - peakspike < tafter  
      % if spike has fewer samples after peak, append NaNs 
      sp = appendfn( sp,     -(nSp - peakspike - tafter) );
   elseif nT - peaktemp < tafter
   end
   
   if all( isnan( sp ) )
      str = sprintf( 'addSpikeToTemplate:error: entire spike is NaN\n' );
      cprintf( 'keywords*', str );
      return;
   end


end