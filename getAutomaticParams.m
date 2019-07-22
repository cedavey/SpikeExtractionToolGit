% GETAUTOMATICPARAMS Returns the ideal parameters for the chosen method
% that fits the input voltage.
%
% Syntax: params = getAutomaticParams(tseries, method);
%
% Inputs:
%     tseries	- The voltage data
%     noisestd - Initial noise standard deviation obtained previously
%     methods	- Method with which the data is being analyzed
%     params   - Default params to keep the structure
%
% Output:
%     params	- Automatic parameters
%
% Artemio - 16/July/2019
function params = getAutomaticParams(tseries, noisestd, method, params)
   str = sprintf('\tStarting automatic parameter configuration...\n');
   printMessage('off','Keywords', str);
   
   % Check noise
   if isempty(noisestd)
      str = sprintf('\tAutomatic parameter configuration failed. Noise input was empty. Using default parameters\n');
      printMessage('off','Errors',str);
      return
   end
   
   % Check data type
   switch lower(tseries.type)
      case 'voltage'
         % Check method
         try
            switch lower(method)
               case 'particle filter'
                  params = calculateAutomaticParams(tseries, noisestd, params);
               case 'recursive least squares'
                  params = calculateAutomaticParams(tseries, noisestd, params);
               otherwise
                  % There is no automatic parameter configuration for the
                  % current data type and method combination.
                  str = sprintf('\tThere is no automatic parameter configuration for the data type: ''%s'' and method: ''%s''\n\tUsing the default parameters\n',tseries.type, method);
                  printMessage('off','Errors',str);
                  return
            end
         catch E
            str = sprintf('\tSomething went wrong with the automatic parameters configuration\n');
            runtimeErrorHandler(E, 'message', str);
            return
         end
      otherwise
         % There is no automatic parameter configuration for the current
         % data type and method combination.
         str = sprintf('\tThere is no automatic parameter configuration for the data type: ''%s'' and method: ''%s''\n\tUsing the default parameters\n',tseries.type, method);
         printMessage('off','Errors',str);
         return
   end
   str = sprintf('\tFinished automatic parameter configuration\n');
   printMessage('off','Keywords', str);
end

% Calculates the parameters for particle filter and recursice least squares
function params = calculateAutomaticParams(tseries, noisestd, params)
   % Check the first 30 spikes in the recording
   spike_amp = []; % Spike amplitude
   spike_smp = []; % Spike sample (time at which the spike happened)
   v = tseries.data;
   dt = tseries.dt;
   th = noisestd * 3;
   idxst = 0;
   idxend = 0;
   smpls = 300;
   peakfn = getPeakFn(params.select_peaks.value);
   
   % Loop until found 30 spikes
   while numel(spike_amp) < 30
      idxst = idxend + 1;
      idxend = idxst + smpls;
      
      warning('off','signal:findpeaks:largeMinPeakHeight');
      [a, b] = findpeaks( peakfn( v(idxst:idxend) ), 'MinPeakHeight', th, 'MinPeakWidth', 3);
      spike_amp = [spike_amp; a];
      spike_smp = [spike_smp; idxst+b-1];
   end
   
   rms = sqrt((sum(spike_amp.^2))/numel(spike_amp));
   varn = noisestd^2; % Variance of the noise
   vars = var(spike_amp); % Variance of the spike amplitudes
   SNR = abs(10*log10(rms/varn));
   sp_rate = numel(spike_amp) / (spike_smp(end) * dt); % Firing rate
   
   % Chose automatic parameters based on spike rate and SNR. The factors
   % are arbitrary, but they can be updated.
   params.voltage_magnitude.value = sqrt(SNR); % Spike threshold is half of the signal to noise ratio
   params.glitch_magnitude.value = 2 * SNR; % Glitch threshold is twice the signal to noise ratio
   params.forgetting_factor.value = 0.99; % This has proved to be the best value
   params.jump_ahead.value = 500/sp_rate;
end

% we're identifying peaks by thresholding, so if we want only positive,
% return actual voltage, if we want negative return the negative of
% voltage, and if we want both positive & negative return abs val.
function peakfn = getPeakFn( select_peaks )
   switch select_peaks
      case 'positive'
         peakfn = @(x) x;
      case 'negative'
         peakfn = @(x) -x;
      case {'all', 'both'}
         peakfn = @(x) abs(x);
   end
end