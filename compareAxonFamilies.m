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

   mutmp1  = fam1.meanspike;   % mean template for family 1
   mutmp2  = fam2.meanspike;   % mean template for family 2
   peak1   = fam1.mup;         % moving average positive peak for family 1
   peak2   = fam2.mup;         %
   trough1 = fam1.mun;         % moving average negative peak for family 1
   trough2 = fam2.mun;
   stdp1   = sqrt(fam1.varp);  % moving average std dev of positive peaks in family 1
   stdp2   = sqrt(fam2.varp);
   stdn1   = sqrt(fam1.varn);  % moving average std dev of negative peaks in family 1
   stdn2   = sqrt(fam2.varn);
   minvar  = min( fam1.varp, fam2.varp );
   alignment_ind1 = fam1.alignment_ind;
   alignment_ind2 = fam2.alignment_ind;
   %% Check if the shapes highly correlated
   % lagged corr btwn template & spike
   [ mutmp1, mutmp2,stimes1, stimes2,alignment_ind1, alignment_ind2  ] = alignDiffLengthSpikes( mutmp1, mutmp2, fam1.stimes, fam2.stimes,alignment_ind1,alignment_ind2, true );
   [ rho, ind ]  = xcorr( mutmp1, mutmp2, 'normalized' );
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

   % I can't use alignDiffLengthSpikes to align the spikes now that they're
   % added to maximise cross correlation rather than aligned at the peak.
   % I need to align at a marker (e.g., the spike peaks), to get rid of
   % potentially huge differences in length. Then do minor alignments to
   % maximise lagged correlation. If the former, we can simply get
   % the indices of the mean, and move all spikes to match this template by
   % aligning their peaks with the template peak. However, this shift will
   % reduce correlation between the spikes, since maximising cross-correlation
   % does not always result in spikes aligned at their peak. Therefore, I am
   % writing 2 functions, one to align at the peak, and the other to simply
   % make the spikes the same length as the template, to ensure that we
   % maintain maximal cross-correlation.

   % The following functions align the family spikes to the templates by
   % aligning their peaks to the template peaks
%    [ spikes1 ]  = alignToMarkerInd( mutmp1, fam1.spikes );
%    [ spikes2 ]  = alignToMarkerInd( mutmp2, fam2.spikes ); 
    % REMOVED BECAUSE THEY ARE ALREADY FULLY ALIGNED WITHIN OWN FAMILY
    
   % The following functions align the family spikes to the templates by
   % shifting to maximise the cross-correlation between the templates
   % Align spikes to maximise the cross-correlation. We only need to shift
   % one family of spikes by index of max lagged corr (i.e., peakind) to 
   % align them. The other family can just be culled to be the same size.
   % Choose the larger family to stay still, & the smaller to shift.
   dt = nanmean( stimes1(2:end,1) - stimes1(1:end-1,1)); % some can be nan so calc. mean
   if fam1.Nspikes > fam2.Nspikes
      spikes2   = alignToMaxCorr( spikes2,  peakind );
      stimes2   = stimes2 - dt*peakind;
       mutmp2   = alignToMaxCorr(  mutmp2,  peakind );
       alignment_ind2 = alignment_ind2+peakind; 
   else
      % gotta take negative of peakind since index was relative to fam1,
      % and tells us how much ahead/behind fam2 is for max corr
      spikes1   = alignToMaxCorr( spikes1, -peakind );
      stimes1   = stimes1 + dt*peakind;
       mutmp1   = alignToMaxCorr(  mutmp1, -peakind );
      alignment_ind1 = alignment_ind1-peakind;
   end
   % Now that we're aligned to max cross correlation, we need to make the
   % spikes in each family a consistent size with the templates


   % Merge the 2 families together now that everything is a consistent size
   Nspikes      = fam1.Nspikes + fam2.Nspikes;
   meanspike    = ( mutmp1 * fam1.Nspikes + mutmp2 * fam2.Nspikes ) / Nspikes; 
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


% Align spikes to a master spike using marker indicies (such as the mean 
% spike's peak, or zero crossing). Make all spikes the same length as the
% given template, plus align them with the peak
% Inputs:
%   inSpikes - unaligned spike matrix of unspecified length, with T X N
%   peakind  - index of max lagged correlation 
function [ spikes ]  = alignToMaxCorr( spikes, peakind )
try % remove this once the code is working
   % Align spikes to maximise the cross-correlation
   [T, N] = size(spikes);   
   if peakind>0
      spikes = [zeros(peakind,N); spikes(1:end-peakind,:)];
   elseif peakind<0
      spikes = [spikes(1-peakind:end,:); zeros(peakind,N)];
   end
   
catch ME
   disp('Problem with the merge');
end
end


% Align spikes to a master spike using marker indicies (such as the mean 
% spike's peak, or zero crossing). Make all spikes the same length as the
% given template, plus align them with the peak
% Inputs:
%   inSpikes - unaligned spike matrix of unspecified length, with T X N
%   tbefore  - number of samples before the marker
%   tafter   - number of samples  after the marker
% (Note: this is all literally copied directly from alignDiffLengthSpikes,
%  with shorten==true)
function [ outSpikes ] = alignToMarkerInd( template, inSpikes )
try % remove this once the code is working
   markerfn    = @getMaxInd;                 % align to peak
   Ttemp       = size(template,1);           % num time samples in template spike
   peakTemp    = markerfn( template );       % index of peak template sample
   tbefore     = peakTemp;                   % samples before marker 
   tafter      = Ttemp - peakTemp; % samples after peak

   [T, N]      = size(inSpikes);             % num time samples x num spikes
   outSpikes   = zeros(Ttemp,N);             % output spikes with same size as template

   % Since we want to max cross-corr, an individual spike may not be
   % aligned with the template peak, so just calculate the mean peak
   % instead, and align the mean peak to the template mean peak (which
   % should be the same, since the template is the mean spike)
   muSpike     = mean( inSpikes, 2 );        % mean family spike
   peakSp      = markerfn(muSpike);          % get index of mean spike's peak 
   sbefore     = peakSp;                     % number of spike samples before peak 
   safter      = T - sbefore;                % number of spike samples after peak

   % the input spikes might be longer or shorter than the template, so
   % grab as many samples as you can before & after the mean spike's peak, 
   % and pop them into the spike vector


   % get samples before peak
   if sbefore > tbefore                      % pre-peak spike is longer than template
      outSpikes(1:peakTemp,:) = inSpikes(1+(sbefore-tbefore):peakSp,:);
   else                                      % spike is shorter than template
      outSpikes(1+(tbefore-sbefore):peakTemp,:) = inSpikes(1:peakSp,:);
   end

   % get samples after peak
   if safter > tafter                     % post-peak spike is longer than template
      outSpikes(peakTemp:end,:) = inSpikes(peakSp:end-(safter-tafter),:); 
   else                                   % post-peak spike is shorter than template
      outSpikes(peakTemp:end-(tafter-safter),:) = inSpikes(peakSp:end,:);
   end

catch ME
   disp('Problem with the merge');
end

end




















