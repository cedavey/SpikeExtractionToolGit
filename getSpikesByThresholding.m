% [spikes, stimes] = getSpikesByThresholding(tseries, params, templateN)
% Extract spikes from time series by thresholding positive & negative
% voltage, requiring the negative peak to be right after the positive peak,
% and imposing a minimum duration on both the positive and negative
% component of the spike. You can provide a required number of samples in
% each spike, if you're going to match them to a template, and the index
% number of the peak. This is so all spikes are aligned else the covariance
% might be low even though they're visually very similar, but out by a tstep
%
% Inputs:
%  tseries     - time series vector
%  params      - params entered via gui
%  templateN   - number of samples in the template spike (optional)
%  peakN       - the peak index of the template spike to allow alignment (optional)
%
% Outputs:
%  spikes      - matrix of voltage x spikes, where we've used the mean spike
%                spike length to extract spikes of all the same length
%  stimes      - cell array of spike times
%  spikesFull  - cell array of spike voltages, with full spike size rather
%                than forcing all spikes to be the same length
%  stimesFull - cell array of spike times, the full spike rather than
%               forcing all to be same length
%  alignment_inds - index in variable length spikes to align on
function [spikes, stimes, spikesFull, stimesFull,alignment_inds] = getSpikesByThresholding(tseries, params, templateN, peakN)
   spikes = []; stimes = []; spikesFull =  []; stimesFull = []; alignment_inds = [];
   dt         = double( tseries.dt );
   voltage    = double( tseries.data(:) );
   time       = double( tseries.time(:) );
   
   str = sprintf( '\tIdentifying possible spikes...\n' );
   printMessage('off', 'Keywords', str );

   % threshold the voltage timeseries - remove values that are too
   % small, and positive values that remain positive for a short period
   posthresh  = params.positive_threshold.value;   % min pos voltage (std dev)
   negthresh  = params.negative_threshold.value;   % min neg voltage (std dev)
   if posthresh > 10
      str = sprintf( '\tWarning: positive threshold for spike extraction is very large (%2.2f)\n', posthresh );
      cprint( 'Keywords*', str );
   end
   if negthresh > 10
      str = sprintf( '\tWarning: negative threshold for spike extraction is very large (%2.2f)\n', negthresh );
      cprint( 'Keywords*', str );
   end
   if isfield( params, 'glitch_threshold') && params.glitch_threshold.value > 0
      glitchthresh = params.glitch_threshold.value;   % glitch threshold
      doglitch= true;
   else
      doglitch= false;
   end
   minpostime = params.min_positive_duration.value;% min duration (ms)
   minnegtime = params.min_negative_duration.value;% min duration (ms)
   avg_window = params.avg_window.value;           % size of avg window to incl in calc (seconds)
   avg_window = round(avg_window) / dt;            % window in seconds, convert to samples
   skip_window= params.skip_window.value;          % how much to skip ahead for next calc (seconds)
   skip_window= round(skip_window) / dt;           % skip window in seconds, convert to samples
   avg_window = ternaryOp(iseven( avg_window),  avg_window,  avg_window+1); % make window even
   skip_window= ternaryOp(iseven(skip_window), skip_window, skip_window+1); % make window even

   minpostime = round((minpostime*1e-3)/dt);       % min duration given in ms
   minnegtime = round((minnegtime*1e-3)/dt);       % min duration given in ms
   
%   % find where voltage goes from -ve to +ve
%   posv       = voltage>0;
%   changePos  = find(diff([0; posv; 0])==1);       % location of change from 0 to 1
%   % find where voltage goes from +ve to -ve
%   changeNeg  = find(diff([0; posv; 0])==-1);      % location of change from 0 to -1
%   if voltage(1)<0
%      timepos = changeNeg - changePos;             % duration of +ve
%      timeneg = [changePos(1)-1; changePos(2:end) - changeNeg(1:end-1)];
%      % spike starts with positive so ditch initial negative
%      timeneg = [timeneg(2:end); 0];
%   else
%      timepos = changeNeg - changePos;          % duration of +ve
%      timeneg = [changePos(2:end) - changeNeg(1:end-1); 0];
%      % since we alternative from pos to neg, at each index position
%      % we have a +ve spike, followed by a -ve spike, so apply thresh
%   end        
%   pass       = timepos>minpostime & timeneg>minnegtime; 
%   if sum(pass)==0
%      if sum( timepos > minpostime ) == 0
%         displayErrorMsg('No spikes have sufficient positive duration, abandoning ship...');
%      end
%      if sum( timeneg>minnegtime ) == 0
%         displayErrorMsg('No spikes have sufficient negative duration, abandoning ship...');
%      end   
%      return;

   % make the template - need to work out how to set this in gui
   getvalence = @(i) params.first_phase_pos.value.*prod(-1.*ones(i+1,1));
   template_params = arrayfun(@(i) struct('voltage_sign',getvalence(i),'min_duration_inds',ternaryOp(getvalence(i)==1,minpostime,minnegtime),...
                                          'max_duration_inds',inf,'minamplitude',ternaryOp(getvalence(i)==1,posthresh,negthresh),...
                                          'do_alignment',i==params.alignment_phase.value),1:params.phasenumber.value);
  
   [Spike_start_ind,Spike_duration,spikes,sinds,stimes,extreme_values,alignment_inds] = threshold_based_on_block_template(voltage,time,template_params);

   if length(spikes) == 0
      displayErrorMsg('No thresholded spikes met positive/negative duration criteria, exiting...');
      return

   end
   [stddev,noise,ntime] = estimate_rolling_noise_std(voltage,time,Spike_start_ind,Spike_duration,avg_window,glitchthresh,params);
   
   %construct a threshold for the blocks in the spike based on the standard
   %deviation of the noise 
   constructhresh = stddev(Spike_start_ind).*horzcat(template_params(:).minamplitude);
   keep = all(extreme_values>=constructhresh,2);
   
   if sum(keep) == 0
      displayErrorMsg('No spikes passed the amplitude criteria., Giving up, but plotting a helpful figure')
      make_helpful_error_plot(template_params,ntime,noise,time,Spike_start_ind,extreme_values,constructhresh,'amplitude')
       return;
   end

   spikes = spikes(keep);
   stimes = stimes(keep);
   sinds  = sinds(keep);
   alignment_inds = alignment_inds(keep);
   % Glitch threshold
   if doglitch
       if avg_window>0
        [ glitch ] = getGlitchThreshold( glitchthresh, extreme_values(keep,:) );
       else
        [ glitch ] = getGlitchThreshold( glitchthresh, extreme_values(keep,:));
       end
      keep_glitch =  all(extreme_values(keep,:) <= glitch,2);%(pos_max(keep) <= glitch_pos) & (neg_min(keep) >= glitch_neg);
      if sum( keep_glitch ) == 0
          %all of the spikes are glitches so none of them are kept
           displayErrorMsg('No spikes survived the positive glitch threshold, sadly giving up but plotting a helpful figure...')
           make_helpful_error_plot(template_params,ntime,noise,time,Spike_start_ind(keep),extreme_values(keep,:),glitch,'glitch') 
           return;
      elseif sum(keep_glitch) < 100
         str = sprintf('Warning: not many spikes (%d) to sort by AP template', sum(keep));
         displayErrorMsg(str);
      end
   
      spikes = spikes(keep_glitch);
      stimes = stimes(keep_glitch);
      sinds  = sinds( keep_glitch);
      alignment_inds = alignment_inds(keep_glitch);
   end
   
   % now make sure spike peaks are aligned else they might be the same but
   % have just slightly different befores/afters - exclude beginning & end
   % part so max isn't at first or last samples rather than middle peak
   offset     = 3;
   
   % we must remove any spikes where the length is less than 2*offset
   remove = cellfun(@length, spikes)<=offset*2;
   spikes(remove) = [];
   stimes(remove) = [];
   sinds(remove) = [];
   alignment_inds(remove) = [];
   

   % user wants template spikes to be a certain length
   if exist('templateN','var') && exist('peakN','var')
      tbefore = peakN;
      tafter  = templateN - peakN; 
      
   % determine best template spike length from avg time before/after peaks
   else
      tbefore = ceil( median(alignment_inds) );
      splength= cellfun( @length, spikes );
      splength= ceil( median(splength) );
      tafter  = ceil( median(splength) - tbefore );
   end
   % get mean time before & after spike peak - extend voltage so don't go over the edge
   voltage    = [voltage; zeros(tafter,1)];
   time       = [ time; time(end)+(1:tafter)'*dt ];
   % if first spike peaks before avg pre-peak time, add zeros to it
   if round( stimes{1}(alignment_inds(1)) / dt ) < tbefore
      voltage = [ zeros(tbefore,1); voltage ];
      time    = [ time(1)-(tbefore:-1:1)' * dt ; time ];
      shift   = tbefore;
   else
      shift   = 0;
   end
   
   spikesFull = spikes;
   stimesFull = stimes;
   sp         = arrayfun( @(i) voltage( sinds{i}(alignment_inds(i))+shift-tbefore+1 : sinds{i}(alignment_inds(i))+shift+tafter ), ...
                              1:length(spikes),'uniformoutput', false );
   spikes     = cell2mat( sp );
   stimes     = arrayfun( @(i) time(    sinds{i}(alignment_inds(i))+shift-tbefore+1 : sinds{i}(alignment_inds(i))+shift+tafter ), ...
                              1:length(sp),'uniformoutput', false );
%    svoltage   = arrayfun( @(i) voltage( sinds{i}(peakind(i))+shift-tbefore+1 : sinds{i}(peakind(i))+shift+tafter ), ...
%                               1:length(spikes),'uniformoutput', false );
                          
                          
   str = sprintf( '\tDone.\n' );
   printMessage('off', 'Keywords', str );

end


function [glitch_thresh, tindex] = getGlitchThreshold(glitchthresh, extreme_values)
   window = 300;
   skip   = 1;
   if size( extreme_values , 1) > window
       
      [avgs,~,tindex] = movingAverageVariance( extreme_values , window, skip, @mean, true);
   else
      tindex = [];
      avgs = nanmean(extreme_values,1);
      
   end
   glitch_thresh = avgs * glitchthresh;    
end

function [stddev,noise,ntime] = estimate_rolling_noise_std(voltage,time,Spike_start_ind,Spike_duration,avg_window,glitchthresh,params)
 
   % assume noise is all the bits that don't have positive & negative
   % consecutive runs - i.e. get from end of last spike to beginning of next
   noise = arrayfun(@(i) toVec( voltage( (Spike_start_ind(i)+Spike_duration(i)+1) : Spike_start_ind(i+1)-1 ) ), ...
                         1:length(Spike_start_ind)-1, 'uniformoutput', false);
   % append noise before first potential spike
   noise = [{toVec( voltage(1:Spike_start_ind(1)) )}, noise(:)'];
   noise = cell2mat( noise(:) );
   % get time of each noise sample - do same as for noise but get values of
   % time instead of voltage signal
   ntime = arrayfun(@(i) toVec( time( (Spike_start_ind(i)+Spike_duration(i)+1) : Spike_start_ind(i+1)-1 ) ), ...
                         1:length(Spike_start_ind)-1, 'uniformoutput', false);
   % append noise before first potential spike
   ntime = [{toVec( time(1:Spike_start_ind(1)) )}, ntime(:)'];
   ntime = cell2mat( ntime(:) );

  
   %can't do a moving average if our sample is too short - so turn it off
   if avg_window >= length(time)
       avg_window = 0;
   end 
   
   % calculate moving avg std dev for thresholding potential spikes
   if avg_window>0 
      
      %There can be glitches in the noise vector remove samples past the
      %glitchthreshold*std(noise) both positive and negative
      noise(noise>glitchthresh*std(noise)|noise<-glitchthresh*std(noise))=0;
       
       
      % only include samples that aren't potential spikes
      [sigma, tsmooth] = movingAverageUnevenSamples( ntime, noise, params.avg_window.value, params.skip_window.value, @std );
      sigma(tsmooth==0) = [];
      tsmooth(tsmooth==0) = [];
      % add first & last samples else interp1 thinks they're 0, which makes
      % beginning & end stddev v.small
      tsmooth = [time(1); tsmooth];
      sigma   = [sigma(1); sigma]; 
      if tsmooth(end) < time(end)
         tsmooth = [tsmooth; time(end)];
         sigma   = [sigma;  sigma(end)];
      end
      stddev  = interp1( tsmooth, sigma, time );
   % no moving avg std dev required, so calculate direction
   else
      stddev  = std(voltage);
   end
   
%    posthresh  = stddev*posthresh;
%    negthresh  = stddev*negthresh; 
% 
%    % get threshold at appropriate point in timeseries if h ave mvg avg std dev
%    if avg_window > 0
%      posthresh= posthresh(startind);
%      negthresh= negthresh(startind);
%    end
end

function [] =  make_helpful_error_plot(template_params,ntime,noise,time,Spike_start_ind,extreme_values,thresh,threshname)
    figure; r = length(template_params); c=1; hold on;
    if size(thresh,2) ==1
        thresh = repmat(thresh,1,r);
    end
      for row = 1:r
          subplot(r,c,row)
          plot( ntime, noise, 'k', 'linewidth', 1 );
          hold on
          plot( time(Spike_start_ind), template_params(row).voltage_sign.*extreme_values(:,row),'.' ); 
          plot( time(Spike_start_ind), template_params(row).voltage_sign.*thresh(:,row), '--' )
          legend( 'noise', 'spike extreme value', [threshname ' threshold'] );
          title(sprintf(' %s of part %d of spike compared to %s threshold',ternaryOp(template_params(row).voltage_sign==1,'maxima','minima'),row,threshname))
          xlabel( 'Time (s)' );
          ylabel( 'Voltage' );
          hold off
       end

end

function [Spike_start_ind,Spike_duration,spikes,sinds,stimes,extreme_values,alignment_inds] = threshold_based_on_block_template(voltage,time,template_params)
%template_params is a struct array detailing each block of the template. Template
%must have at least two blocks (each represented by element in struct
%array, with following information
% voltage_sign: 1 or -1 (adjacent blocks must have alternating voltage
% signs
% min_duration_inds: integer - minimum number of time bins spent
% postive/negative
% max_duration_inds: integer - maximum number of time bins spent
% positive/negative
% minamplitude: noisestd minimum to pass ( always positive because will
% multiply by voltage_sign
   posv = voltage>0;
   Nblocks = length(template_params);
   Startblocksign = template_params(1).voltage_sign;
   Endblocksign = template_params(end).voltage_sign;
   
   zerocross = find(abs(diff(posv))>0);
   
   if (Startblocksign == 1) == posv(1)
       %if the first segment of the time series is the valence we're
       %looking for include 1 as a zerocrossing
       zerocross = [1;zerocross];
   end
   
   if (Endblocksign ==1 ) == posv(end)
       %if the final segment of recording matches what we're looking
       %for, add the length of the timeseries as a zero crossing
       zerocross = [zerocross;length(posv)];
   end
   
   % get durations of each block
   dur = diff(zerocross)';
   block_durations = arrayfun(@(b) {dur(b:2:end-Nblocks+b)},1:Nblocks);
   block_durations = vertcat(block_durations{:});
   start_inds = zerocross(1:2:end-Nblocks+1)+1;
      
   %only find spikes that pass the duration tests for all blocks
   pass = all((block_durations>vertcat(template_params(:).min_duration_inds)) & (block_durations<vertcat(template_params(:).max_duration_inds)),1);
   
   if sum(pass) == 0
       %no spikes found, exit
      [Spike_start_ind,Spike_duration,spikes,sinds,stimes,extreme_values] = deal([]);
      return;
   end
   
   Spike_start_ind = start_inds(pass==1);
   Spike_duration = sum(block_durations(:,pass==1),1);
   Spike_block_durations = block_durations(:,pass==1);
   
   
   
   %Now for each spike, we want to provide an extreme value for each block
   %(minimum if its negative, maximum if it's positive
   spikes     = arrayfun(@(i) toVec(voltage(Spike_start_ind(i):(Spike_start_ind(i)+Spike_duration(i)-1))), ...
                              1:length(Spike_start_ind), 'uniformoutput', false); 
   sinds      = arrayfun(@(i) toVec( Spike_start_ind(i):(Spike_start_ind(i)+Spike_duration(i)-1)), ...
                              1:length(Spike_start_ind), 'uniformoutput', false); 
   stimes     = arrayfun(@(i) toVec(time(Spike_start_ind(i):(Spike_start_ind(i)+Spike_duration(i)-1))), ...
                              1:length(Spike_start_ind), 'uniformoutput', false); 
   
                          
   block_inds = @(b,s) 1+sum(Spike_block_durations(1:b-1,s)):sum(Spike_block_durations(1:b,s));
   extreme_values  = arrayfun(@(s) {arrayfun(@(b) nanmax(template_params(b).voltage_sign.*spikes{s}(block_inds(b,s))),1:Nblocks)} ,1:length(spikes));
   extreme_values = vertcat(extreme_values{:});
   
   % Find the index of the peak to align on
   ab = find([template_params(:).do_alignment],1);
   get_alignment_ind =  @(s,bl,offset,ab) sum(bl(1:ab-1)) + getMaxInd(template_params(ab).voltage_sign.*s(sum(bl(1:ab-1))+1:sum(bl(1:ab))));
   
   alignment_inds = arrayfun(@(s) get_alignment_ind(spikes{s},Spike_block_durations(:,s),0,ab) ,1:length(spikes));
   
   
   
end