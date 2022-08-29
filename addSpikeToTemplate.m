% [ temp, tempAPs ] = addSpikeToTemplate( temp, sp, tempAPs ) 
%
% [ temp, tempAPs, sptime, tempTimes ] = addSpikeToTemplate( temp, sp, tempAPs, sptime, tempTimes )
%
% for spikes and templates that allow different lengths, when you add a new
% spike to a template you need to align all the peaks, and then take an
% ensemble average where the number of samples contributing to each time
% point can be different, depending on the lengths of the constituent spikes
% Inputs:
%  temp     - template vector containing a single AP shape
%  sp       - spike vector to be merged with existing template
%  tempAPs  - matrix of spikes contributing to the template, aligned at the
%             peak & with NaNs at beginning or end if a particular spike
%             had no samples at this point
%  sptimes  - vector of spike times
%  tempTimes - time vector for spike templates
function [ temp, tempAPs, sptime, tempTimes ] = addSpikeToTemplate( temp, sp, tempAPs, sptime, tempTimes )
   if nargin==5, doTime = true; else, doTime = false; end
   peakind  = getMaxInd( temp ); % find index of template peak 
   peaksp   = getMaxInd( sp   ); % find index of spike peak 

   nSp      = length( sp   );    % length of spike
   nT       = length( temp );    % length of template
   tbefore  = max( peakind, peaksp );
   tafter   = max( nT - peakind, nSp - peaksp );
   
   % append NaNs to template spikes or new spike to make them the same size
   prependfn= @( sp, N ) [ NaN( N, size(sp,2)); sp ];
   appendfn = @( sp, N ) [ sp; NaN( N, size(sp,2)) ];
   if peaksp < tbefore 
      % if spike has fewer samples before peak, prepend NaNs 
      sp = prependfn( sp, tbefore - peaksp );
      if doTime
         sptime = prependfn( sptime, tbefore - peaksp );
      end
   elseif peakind < tbefore
      % if template constituent spikes have fewer samples before peak, append 
      % NaNs to all of the constituent spikes 
      tempAPs = prependfn( tempAPs,   tbefore - peakind );
      if doTime
         tempTimes = prependfn( tempTimes, tbefore - peakind );
      end
   end
   if nSp - peaksp < tafter  
      % if spike has fewer samples after peak, append NaNs 
      sp = appendfn( sp,     -(nSp - peaksp - tafter) );
      if doTime
         sptime = appendfn( sptime, -(nSp - peaksp - tafter) );
      end
   elseif nT - peakind < tafter
      % if template constituent spikes have fewer samples before peak, append 
      % NaNs to all of the constituent spikes 
      tempAPs   = appendfn( tempAPs,   -(nT - peakind - tafter) );
      if doTime
         tempTimes = appendfn( tempTimes, -(nT - peakind - tafter) );
      end
   end
   
   if all( isnan( sp ) )
      str = sprintf( 'addSpikeToTemplate:error: entire spike is NaN\n' );
      cprintf( 'keywords*', str );
      return;
   end
   % now add the new spike to the constituent spikes for the template
   try
      tempAPs   = [tempAPs sp];
      if doTime
         tempTimes = [tempTimes sptime];
      end
   catch ME
      str = sprintf( 'addSpikeToTemplate:error: Error concatenating spikes\n' );
      cprintf( 'keywords*', str );
      runtimeErrorHandler( ME );
   end
   % now take an ensemble average of template, where the number of samples
   % contributing to each time point may change across the spike, depending
   % on the length of the constituent spikes 
   % Note: when taking the average, while we keep the entire spike for all
   % of the constituent spikes, the templates are getting fucked up by
   % increasingly wonky start and ends because if one spike has a wacky
   % extra beginning bit, then the ensemble average there is just from it,
   % so the template goes a bit wacky, so then it matches the next slightly
   % wackier spike, and on it goes, until the whole template is wacky
   perc  = 2/3; % need at least this percent of spikes to have a sample at each pt
   temp  = nanmean( tempAPs, 2 );
%    [temp, tempAPs] = cleanTemplates( temp, tempAPs, perc );
   
   nAPs  = size( tempAPs, 2 );
   Nreqd = ceil( nAPs * perc );
   if nAPs > Nreqd
      % only do it randomly so there's a chance to build up some spikes
      if rand(1) > 0.8
         Nensemble = sum( ~isnan( tempAPs ), 2 );
         ditch = Nensemble <= Nreqd;
         if sum( ditch ) == size( temp, 1 )
            str = sprintf( 'Attempting to ditch %d non-overlapping samples of %d template samples\n', ...
                            sum( ditch ), sum( ~ditch ) );
            cprintf( 'Keyword*', str );
         else
            temp( ditch ) = [];
            tempAPs( ditch, : )   = [];
            if doTime
               tempTimes( ditch, : ) = [];
            end
         end
      end
   end

   if isempty( temp ) || isempty( tempAPs )
      str = sprintf( 'Adding spike to template has decimated it...moron\n' );
      cprintf( 'Keyword*', str );
   end
end
