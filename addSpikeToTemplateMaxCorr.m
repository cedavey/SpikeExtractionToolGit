% Add current spike to template by making the template, family spikes, and
% current spike, all the same length, and then optimising lagged
% correlation. Spike times for the family spikes and current spike are
% updated to have the same length vectors also. 
function [ temp, tempAPs, sp, sptime, tempTimes,alignment_ind_template,alignment_ind_spike] = addSpikeToTemplateMaxCorr( ...
                                                         template, curr_spike,alignment_ind_template,alignment_ind_spike, tempAPs, sptime, tempTimes, lambda )
   % We can either cull the longer spike, or append 0s to the shorter spike
   % (Note that we treat the samples before and after the peak or zero
   % crossing separately, so there will be 2 culling/appending events)
   shorten = false;       % shorten spikes to equal length, or lengthen spikes to equal length

   % align at peak first to make sure we're not using wiggly beginning or end bits
   tempcopy = template; spcopy = curr_spike; timecopy = tempTimes; sptimecopy = sptime; tempAPcopy = tempAPs;
%    [ template, curr_spike, tempTimes, sptime, tempAPs,pad_inds] = alignDiffLengthSpikes( template, curr_spike, tempTimes, sptime, shorten, tempAPs ,(1:length(template))');
%    pad_inds = pad_inds==0;
%    if size(template,1) ~= size(curr_spike,1) || size(template,1)   ~= size(tempAPs,1) ...
%       size(sptime,1)    ~= size(template,1)  || size(curr_spike,1) ~= size(tempTimes,1)
%       disp('problem with sizes');
%    end
   % anchor to max corr when using lagged corr, rather than spike peaks
   nT      = numel(template); 
   nS      = numel(curr_spike); 
   nAP     = size(tempAPs,2);
   
   % extra time added to template of spike needs to be accounted for in
   % time vector
   dt      = nanmean(sptime(2:end) - sptime(1:end-1));
   addtime = (1:(abs(nT-nS)) )' *dt;

   temp    = template;
   sp      = curr_spike; 
   % pad with zeros if one is smaller than the other - if a lagged corr is
   % the optimal, make sure to move in the direction of overwriting 0s...
   % This correction is meaningless since we've already aligned the spikes
   % earlier - if remove that alignement, will this solve the problem
   if nT > nS
      sp         = zeros(nT,1); 
      sp(1:nS)   = curr_spike; 
      temp       = template; 
      sptime     = [sptime; sptime(end) + addtime];
   elseif nS > nT
      temp       = zeros(nS,1); 
      temp(1:nT) = template; 
      % gotta extend the template also
      allspikes  = zeros(nS,nAP);
      allspikes(1:nT,:) = tempAPs;
      tempAPs    = allspikes; 
      sp         = curr_spike; 
      tempTimes  = [tempTimes; tempTimes(end,:) + addtime];
   end

   % lagged corr btwn template & spike
   [ rho, ind ]  = xcorr( temp, sp, 'normalized' );
   [ rho, peak ] = max(rho);
   % if lag is positive, move spike forward, if negative, move temp forward
   peakind      = ind(peak);
    
try
   % Align spikes to maximise the cross-correlation. 
   
   %are we correctly adjustign the time functions?
   if peakind>0
        sp      = [zeros(peakind,1);   sp(1:end-peakind)];
        sptime =  sptime-peakind*dt;
        alignment_ind_spike = alignment_ind_spike+peakind;
   elseif peakind<0
        temp    = [zeros(-peakind,1); temp(1:end+peakind)];
        tempAPs = [zeros(-peakind,size(tempAPs,2)); tempAPs(1:end+peakind,:)];
        tempTimes = tempTimes+peakind*dt;
        alignment_ind_template = alignment_ind_template - peakind;
   end

   tempAPs    = [tempAPs sp]; % add spike to template's spikes
   Nspikes    = size(tempAPs,2);
   % efficient calculation of mean template
   temp       = recursiveUpdate( temp, [], sp, lambda );
%    template  = ( template * (Nspikes-1) + curr_spike ) / Nspikes;
%    template  = nanmean( tempAPs, 2 );
   tempTimes  = [tempTimes sptime];
%FInd peak nearby old peak location (incase the update has shifted
%it a few indices away
   alignment_ind_template = getMaxInd(temp*sign(temp(alignment_ind_template)).*(((1:length(temp))'-alignment_ind_template)<5));

catch ME
   disp('problem');
end

end


