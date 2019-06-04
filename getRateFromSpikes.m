% [rates, time] = getRateFromSpikes(tseries, method, method_params)
function [rates, time, dt] = getRateFromSpikes(tseries, method, method_params)
   % for each AP template we have a group of families of various scales, in
   % terms of amplitude - get rate for each family within each template
   nTemplates  = length( tseries.data );
   rates       = cell( nTemplates, 1 );
   T           = tseries.time(end);
   dt          = tseries.dt;
   nT          = round( T/dt ); % number of time samples
   
   avg_window  = method_params.moving_avg_window.value;
   skip_window = method_params.skip_window.value;
   avg_samples = round( avg_window / dt );
   skip_window = round( skip_window / dt );
   % convolve spikes with kernel to smooth rate estimate
   kernel_type = method_params.smoothing_kernel.value; 
   nS          = round( T/dt/skip_window ); % number of rate estimates
   if nS < 5 % require at least 5 rate estimates to go ahead
      displayErrorMsg('Oops, the window size to average spike rate over is too long for the time series...');
      rates = []; time = []; dt = []; 
      return;
   end
   
   sig  = 3;
   nflt = max( min( nS, sig*5 ), 3);
   x    = (-nflt:nflt)';
   filt = true;
   switch kernel_type
      case 'boxcar'
         kernel = ones( size(x) );
      case 'exponential'
         kernel = exp( -x / sig );
      case 'gaussian'
         kernel = normpdf(x, 0, sig);

      case 'none'
         kernel = 1;
         filt   = false;
         
   end
   kernel = kernel./norm(kernel);

   % gotta decide whether to avg spikes & then convolve, or convolve with
   % raw spike train
   
   % for each family within the templates we need to find rate - the spike
   % times are given for each family within a cell array
   APstimes = tseries.APstimes;
   
   % hack! when I make timeseries out of consecutive spikes, it means
   % the spiketimes are much later than the length of the timeseries -
   % check that this isn't happening and fucking us up
   latest = ( max( cellfun(@(t) max(cellfun(@(d) max(d(:)), t)), tseries.APstimes) ) / dt );
   if latest > nT
      nT = latest;  % number of samples
      T  = nT * dt; % end time
      nS = round( T/dt/skip_window ); % number of rate estimates
   end
   
   for ti=1:nTemplates
      % Now we have a list of spike times for each AP template, get spike rate
      % by averaging the number of spikes in each window over time. To ensure
      % spike rate is causal, calculate rate at current time point, by
      % averaging from current time to number of seconds prior to current time
      nFam      = length( tseries.APstimes{ti} );
      rates{ti} = zeros( nS, nFam );
      tmp_rates = zeros( nS, nFam );
      
      for ff=1:nFam
         spiketimes = round( APstimes{ti}{ff} / dt );
         
         downt      = 1; % downsampled time index
         for eind = skip_window:skip_window:nT
            sind    = max([eind-avg_samples, 1]);
            nspikes = sum( sind<spiketimes & spiketimes<=eind ); 
            % count spikes over averaging period, & convert to rate
            tmp_rates(downt,ff) = nspikes / avg_window; 
            downt   = downt + 1;
         end
      end
      
      if filt
         smth_rates = zeros( size(tmp_rates) );
         o = conv( ones(size(tmp_rates,1),1), kernel, 'same' );
         for ff=1:nFam
           smth_rates(:,ff) = conv( tmp_rates(:,ff), kernel, 'same' ) ./ o;
         end
         rates{ti} = smth_rates;
         
      else
         rates{ti} = tmp_rates;
      end
   end
   
   % calculate down sampled time vector
   time = (skip_window:skip_window:nT) * dt;
   dt   = skip_window * dt;
end
   
   
   
