% [ sp1, sp2 ] = alignDiffLengthSpikes( sp1, sp2 )
%
% for spikes that are different lengths, this function aligns them at the
% peak and then gets the min number of samples of either spikes before the 
% peak, and the min number of samples of either spike after the peak, to 
% ensure both spikes have valid samples when aligned
function [ sp1, sp2 ] = alignDiffLengthSpikes( sp1, sp2 )
   % if spikes & templates can have diff lengths, use the smaller length
   % to calculate the match, but centre both around the peak 
   nSp     = length( sp1 );
   nF      = length( sp2 ); % length of each template
   peakfsp = getMaxInd( sp2 );
   peaksp  = getMaxInd( sp1 );
   tbefore = min( peakfsp, peaksp );
   tafter  = min( nF - peakfsp, nSp - peaksp );
   
   sp1 = sp1(  peaksp-tbefore+1 :  peaksp+tafter );
   sp2 = sp2( peakfsp-tbefore+1 : peakfsp+tafter );
end