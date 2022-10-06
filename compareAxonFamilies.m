% Determine if 2 axon families have very similar sizes & shapes, and merge
% them if they do

function [similar, mergedfam] = compareAxonFamilies( fam1, fam2, rhothresh, kappa_pos, kappa_neg )
   % fam.meanspike
   % fam.spikes
   % fam.stimes 
   % fam.mup
   % fam.mun
   % fam.varp
   % fam.varn

   % initialise output as failed comparison
   similar   = false; 
   mergedfam = [];

   temp1   = fam1.meanspike; % template for family 1
   temp2   = fam2.meanspike; 
   peak1   = fam1.mup;
   peak2   = fam2.mup; 
   trough1 = fam1.mun; 
   trough2 = fam2.mun;
   stdp1   = sqrt(fam1.varp);
   stdp2   = sqrt(fam2.varp);
   stdn1   = sqrt(fam1.varn);
   stdn2   = sqrt(fam2.varn);
   minvar  = min( fam1.varp, fam2.varp );

   %% Check if the shapes highly correlated
   % lagged corr btwn template & spike
   [ temp1, temp2, stimes1, stimes2 ] = alignDiffLengthSpikes( temp1, temp2, fam1.stimes, fam2.stimes, true );
   [ rho, ind ]  = xcorr( temp1, temp2, 'normalized' );
   [ rho, peak ] = max(rho);
   % if lag is positive, move spike forward, if negative, move temp forward
   peakind      = ind(peak);
   
   if rho <= rhothresh
      return;
   end

   %% Check if the magnitudes are very similar
%    peaks    = abs( peak1   - peak2 )   / max( peak1, peak2 );
%    troughs  = abs( trough1 - trough2 ) / max( trough1, trough2 );
%    compare  = @(s1, s2) abs( max(s1) - max(s2) ) + abs( min(s1) - min(s2) ); % wgt peak & trough the same
%    compare( fam1.meanspike, fam2.meanspike );


   valid1 = ( trough2 >= ( trough1 - kappa_neg * stdn1 ) ) ...
         && ( trough2 <= ( trough1 + kappa_neg * stdn1 ) );
   valid2 = ( peak2   >= (   peak1 - kappa_pos * stdp1 ) ) ...
         && ( peak2   <= (   peak1 + kappa_pos * stdp1 ) ); 

   valid3 = ( trough1 >= ( trough2 - kappa_neg * stdn2 ) ) ...
         && ( trough1 <= ( trough2 + kappa_neg * stdn2 ) );
   valid4 = ( peak1   >= (   peak2 - kappa_pos * stdp2 ) ) ...
         && ( peak1   <= (   peak2 + kappa_pos * stdp2 ) ); 

   % if not similar peak & trough (i.e., within number of std devs), then exit
   if ~( valid1 && valid2 && valid3 && valid4 )
      return;
   end

   %% Merge similar families
   % If we're here then the families have similar size & shape, so merge

   % make all spikes the same length
   [ spikes1, spikes2 ] = alignDiffLengthSpikes( fam1.spikes, fam2.spikes, [], [], true );

   Nspikes      = fam1.Nspikes + fam2.Nspikes;
   meanspike    = ( temp1 * fam1.Nspikes + temp2 * fam2.Nspikes ) / Nspikes; 
   ind          = getMaxInd(meanspike); % get index of max peak for sorting spikes
   peaktimes1   = stimes1(ind,:);
   peaktimes2   = stimes2(ind,:);
   alltimes     = [stimes1 stimes2];
   allspikes    = [spikes1 spikes2];

   [~,sortind]  = sort( [peaktimes1 peaktimes2] );
   stimes       = alltimes(:, sortind);
   spikes       = allspikes(:,sortind);
   curr_spike   = spikes(:,end);
   curr_stime   = stimes(:,end);
   % don't account for number of spikes in each family for the means, since
   % they're recurively updated so older spikes don't contribute much
   mun          = ( fam1.mun  + fam2.mun ) / 2;
   mup          = ( fam1.mup  + fam2.mup ) / 2;
   varn         = ( fam1.varn + fam2.varn ) / 2;
   varp         = ( fam1.varp + fam2.varp ) / 2;

   mergedfam.spikes    = spikes;
   mergedfam.stimes    = stimes;
   mergedfam.meanspike = meanspike;
   mergedfam.last_peak = [min(curr_spike) max(curr_spike)]; 
   mergedfam.last_time = curr_stime( ind ); % time of last spike
   mergedfam.mun       = mun; 
   mergedfam.mup       = mup; 
   mergedfam.varn      = varn;
   mergedfam.varp      = varp;
   mergedfam.Nspikes   = Nspikes;

   % if debugging, we'll have vectors for the recursively updated vars
   if isfield(fam1, 'mup_vector') && isfield(fam2, 'mup_vector')
      mup_vector  = [fam1.mup_vector fam2.mup_vector];
      mup_vector  = mup_vector(sortind);
      mun_vector  = [fam1.mun_vector fam2.mun_vector];
      mun_vector  = mun_vector(sortind);
      varp_vector = [fam1.varp_vector fam2.varp_vector];
      varp_vector = varp_vector(sortind);
      varn_vector = [fam1.varn_vector fam2.varn_vector];
      varn_vector = varn_vector(sortind);

      mergedfam.mup_vector  = mup_vector;
      mergedfam.mun_vector  = mun_vector;
      mergedfam.varp_vector = varp_vector;
      mergedfam.varn_vector = varn_vector;
   end

   similar = true;

   if size(mergedfam.spikes,1) ~= size(mergedfam.stimes,1) ...
   || size(mergedfam.spikes,1) ~= size(mergedfam.meanspike,1)
      disp('Problem with sizes after family merge');
   end
end









