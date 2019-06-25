% [spikes, stimes] = getSpikesByThresholding(tseries, params, templateN)
% Extract spikes from time series by thresholding positive & negative
% voltage, requiring the negative peak to be right after the positive peak,
% and imposing a minimum duration on both the positive and negative
% component of the spike. You can provide a required number of samples in
% each spike, if you're going to match them to a template, and the index
% number of the peak. This is so all spikes are aligned else the covariance
% might be low even though they're visually very similar, but out by a tstep
% Inputs:
%  tseries   - time series vector
%  params    - params entered via gui
%  templateN - number of samples in the template spike (optional)
%  peakN     - the peak index of the template spike to allow alignment (optional)
function [spikes, stimes, sindices] = getSpikesByThresholding(tseries, params, templateN, peakN)
   spikes = []; stimes = [];
   dt         = double( tseries.dt );
   voltage    = double( tseries.data(:) );
   time       = double( tseries.time(:) );
   
   str = sprintf( '\tIdentifying possible spikes...\n' );
   printMessage('off', 'Keywords', str );

   % threshold the voltage timeseries - remove values that are too
   % small, and positive values that remain positive for a short period
   posthresh  = params.positive_threshold.value;   % min pos voltage (std dev)
   negthresh  = params.negative_threshold.value;   % min neg voltage (std dev)
   minpostime = params.min_positive_duration.value;% min duration (ms)
   minnegtime = params.min_negative_duration.value;% min duration (ms)
   avg_window = params.avg_window.value;           % size of avg window to incl in calc (seconds)
   avg_window = round(avg_window) / dt;            % window in seconds, convert to samples
   skip_window= params.skip_window.value;          % how much to skip ahead for next calc (seconds)
   skip_window= round(skip_window) / dt;           % skip window in seconds, convert to samples
   avg_window = ternaryOp(iseven(avg_window), avg_window, avg_window+1); % make window even
   skip_window= ternaryOp(iseven(skip_window), skip_window, skip_window+1); % make window even

   minpostime = round((minpostime*1e-3)/dt);       % min duration given in ms
   minnegtime = round((minnegtime*1e-3)/dt);       % min duration given in ms
   
   % find where voltage goes from -ve to +ve
   posv       = voltage>0;
   changePos  = find(diff([0; posv; 0])==1);    % location of change from 0 to 1
   % find where voltage goes from +ve to -ve
   changeNeg  = find(diff([0; posv; 0])==-1);   % location of change from 0 to 1
   if voltage(1)<0
      timepos = changeNeg - changePos;          % duration of +ve
      timeneg = [changePos(1)-1; changePos(2:end) - changeNeg(1:end-1)];
      % spike starts with positive so ditch initial negative
      timeneg = [timeneg(2:end); 0];
   else
      timepos = changeNeg - changePos;          % duration of +ve
      timeneg = [changePos(2:end) - changeNeg(1:end-1); 0];
      % since we alternative from pos to neg, at each index position
      % we have a +ve spike, followed by a -ve spike, so apply thresh
   end        
   pass       = timepos>minpostime & timeneg>minnegtime; 
   if sum(pass)==0
      displayErrorMsg('No spikes have suffic positive and negative duration, abandoning ship...');
      spikes = []; stimes = []; sindices = []; % This line to prevent run time error "not assigned during getSpikesByThresholding"
      return;
   end
   changePos  = changePos(pass);
   duration   = timepos(pass) + timeneg(pass);
   
   % assume noise is all the bits that don't have positive & negative
   % consecutive runs - i.e. get from end of last spike to beginning of next
   noise = arrayfun(@(i) toVec( voltage( (changePos(i)+duration(i)+1) : changePos(i+1)-1 ) ), ...
                         1:length(changePos)-1, 'uniformoutput', false);
   % append noise before first potential spike
   noise = [{toVec( voltage(1:changePos(1)) )}, noise(:)'];
   noise = cell2mat( noise(:) );
   % get time of each noise sample - do same as for noise but get values of
   % time instead of voltage signal
   ntime = arrayfun(@(i) toVec( time( (changePos(i)+duration(i)+1) : changePos(i+1)-1 ) ), ...
                         1:length(changePos)-1, 'uniformoutput', false);
   % append noise before first potential spike
   ntime = [{toVec( time(1:changePos(1)) )}, ntime(:)'];
   ntime = cell2mat( ntime(:) );

   % now make sure that positive & negative periods are consecutive -
   % get sample right after last negative so we get the final upstroke
   voltage    = [voltage; 0]; time = [time; time(end)+dt]; % so extra sample doesn't break
   spikes     = arrayfun(@(i) toVec(voltage(changePos(i):(changePos(i)+duration(i))+1)), ...
                              1:length(changePos), 'uniformoutput', false); 
   % now make sure that positive & negative periods are consecutive
   sinds      = arrayfun(@(i) toVec( changePos(i):(changePos(i)+duration(i)+1)), ...
                              1:length(changePos), 'uniformoutput', false); 
   stimes     = arrayfun(@(i) toVec(time(changePos(i):(changePos(i)+duration(i)+1))), ...
                              1:length(changePos), 'uniformoutput', false); 
   pos_max    = toVec(cellfun(@nanmax, spikes));
   neg_min    = toVec(cellfun(@nanmin, spikes));
   startind   = toVec(cellfun(@(st) round(st(1)/dt), stimes));
   time       = time(1:end-1); % remove extra sample of time added above
   
   % calculate moving avg std dev for thresholding potential spikes
   if avg_window>0
      % only include samples that aren't potential spikes
      [sigma, tsmooth] = movingAverageUnevenSamples( ntime, noise, params.avg_window.value, params.skip_window.value, @std );
      % add first & last samples else interp1 thinks they're 0, which makes
      % beginning & end stddev v.small
      tsmooth = [time(1); tsmooth];
      sigma   = [sigma(1); sigma]; 
      if tsmooth(end) < time(end)
         tsmooth = [tsmooth; time(end)];
         sigma   = [sigma;  sigma(end)];
      end
      stddev  = interp1( tsmooth, sigma, time );
%       % interp1 does not like having to extrapolate if times of smaller &
%       % larger vectors do not match at beginning and end
%       if sum(isnan(sigma))>0
%          sigma = interp1( tsmooth/(tsmooth(end)-tsmooth(1))*(time(end)-time(1)), sigma, time );
%          sigma(isnan(sigma)) = nanmean(sigma);
%       end

   % no moving avg std dev required, so calculate direction
   else
      stddev  = std(voltage);
   end
   posthresh  = stddev*posthresh;
   negthresh  = stddev*negthresh;

   % get threshold at appropriate point in timeseries if have mvg avg std dev
   if avg_window > 0
     posthresh= posthresh(startind);
     negthresh= negthresh(startind);
   end
   keep       = (pos_max >= posthresh) & (neg_min <= -negthresh);
   if sum(keep)==0
      displayErrorMsg('No spikes have suffic positive and negative amplitude, sadly giving up...');
      return;
   elseif sum(keep)<100
      str = sprintf('Warning: not many spikes (%d) to sort by AP template', sum(keep));
      displayErrorMsg(str);
   end
   spikes     = spikes(keep);
   stimes     = stimes(keep);
   sinds      = sinds(keep);
   
   % now make sure spike peaks are aligned else they might be the same but
   % have just slightly different befores/afters - exclude beginning & end
   % part so max isn't at first or last samples rather than middle peak
   offset     = 3;
   peakind    = getMaxInd( spikes, 1, offset );
   % user wants template spikes to be a certain length
   if exist('templateN','var') && exist('peakN','var')
      tbefore = peakN;
      tafter  = templateN - peakN; 
      
   % determine best template spike length from avg time before/after peaks
   else
      tbefore = ceil(mean(peakind));
      splength= cellfun(@length, spikes);
      splength= ceil(mean(splength));
      tafter  = ceil(mean(splength) - tbefore);
   end
   % get mean time before & after spike peak - extend voltage so don't go over the edge
   voltage    = [voltage; zeros(tafter,1)];
   time       = [time; time(end)+(1:tafter)'*dt];
   % if first spike peaks before avg pre-peak time, add zeros to it
   if round(stimes{1}(peakind(1))/dt) < tbefore
      voltage = [zeros(tbefore,1); voltage];
      time    = [time(1)-(tbefore:-1:1)'*dt; time];
      shift   = tbefore;
   else
      shift   = 0;
   end
   sp         = arrayfun(@(i) voltage( sinds{i}(peakind(i))+shift-tbefore+1 : sinds{i}(peakind(i))+shift+tafter), ...
                              1:length(spikes),'uniformoutput', false);
   spikes     = cell2mat(sp);
   stimes     = arrayfun(@(i) time( sinds{i}(peakind(i))+shift-tbefore+1 : sinds{i}(peakind(i))+shift+tafter), ...
                              1:length(sp),'uniformoutput', false);
   sindices   = arrayfun(@(i) voltage( sinds{i}(peakind(i))+shift-tbefore+1 : sinds{i}(peakind(i))+shift+tafter), ...
                              1:length(spikes),'uniformoutput', false);
                          
                          
   str = sprintf( '\tDone.\n' );
   printMessage('off', 'Keywords', str );
end





