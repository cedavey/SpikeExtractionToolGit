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
   Tbefore  = peaktemp-1;        % time before template peak
   Tafter   = nT - peaktemp;     % time after template peak
   Sbefore  = peakspike-1;       % time before spike peak 
   Safter   = nSp - peakspike;   % time after spike peak

   % append NaNs to template spikes or new spike to make them the same size
   prependfn= @( sp, N ) [ NaN( N, size(sp,2)); sp ];
   appendfn = @( sp, N ) [ sp; NaN( N, size(sp,2)) ];

   % we're modifying the spike to match the template, not the other way
   % around

   % spike peak occurs before template peak
   if Sbefore < Tbefore 
      % if spike has fewer samples before peak, prepend NaNs 
      sp = prependfn( sp, Tbefore - Sbefore );
   elseif Sbefore > Tbefore
      sp(1:(peakspike - peaktemp)) = [];
   end
   if Safter < Tafter  
      % if spike has fewer samples after peak, append NaNs 
      sp = appendfn( sp, -(nSp - peakspike - Tafter) );
   elseif Safter > Tafter
      sp(end-(Safter-Tafter)+1:end) = [];
   end
   
   if all( isnan( sp ) )
      str = sprintf( 'addSpikeToTemplate:error: entire spike is NaN\n' );
      cprintf( 'keywords*', str );
      return;
   end


end