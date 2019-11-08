% voltage = denoiseVoltage(tseries, method, method_params)
% 
% Remove noise from voltage using user specified method from gui
function denoised_voltage = denoiseVoltage(tseries, method, params)
   if nargin==0, help denoiseVoltage; return; end
   denoised_voltage = [];
   % Print denoising message
   str = sprintf( '\tDenoising voltage using %s...\n', method);
   cprintf( 'Keywords', str );
   switch lower(method)
      case 'filter'
%          try
            N   = params.filter_order.value; % number of filter coeffs to use
            hc  = params.high_cutoff_frequency.value; % upper freq bound
            lc  = params.low_cutoff_frequency.value;  % lower freq bound
            df  = 1/tseries.dt;
            
            % Filter functions expect normalised frequency (Hz, not
            % radians) so divide freq by max freq & multiply by 2, so that
            % 1==Nyquist frequency.
            % if low pass freq is 0, design a lowpass filter instead
            if lc==0
               Wc   = hc / df * 2; % cut-off normalised so Nyquist is 1
               type = 'low';
               Wl   = 0; Wh = Wc;  % just in case matlab returns NaNs & gotta design ourself
            elseif hc==0
               Wc   = lc / df * 2; % cut-off normalised so Nyquist is 1
               type = 'high';
               Wl   = Wc; Wh = 0;  % just in case matlab returns NaNs & gotta design ourself
            else
               Wc   = [lc/df*2  hc/df*2];  % cut-off normalised so Nyquist is 1
               type = 'bandpass';
               Wl   = Wc(1); Wh = Wc(2);  % just in case matlab returns NaNs & gotta design ourself
            end
            [b, a]  = butter(N, Wc, type);
            
            % sometimes matlab filter fns fuck up & return NaNs so have to
            % design the filter more manually
            if any(isnan(b)) || any(isnan(a))
               filt = design(fdesign.bandpass('N,F3dB1,F3dB2', N, Wl, Wh), 'butter');
               denoised_voltage = filter(filt, tseries.data);
            
            % this filtering is done in c so is faster
            else
               denoised_voltage = FiltFiltM(b, a, tseries.data);
            end
%          catch ME
%             displayErrorMsg('You''ve fucked something up - check denoiseVoltage:filter');
%          end
         
      case 'threshold'
         try
            [spikes, stimes] = getSpikesByThresholding(tseries, params);
            denoised_voltage = zeros(size(tseries.voltage)); 
            sindices         = round(cell2mat(stimes(:))/dt);
            denoised_voltage(sindices) = cell2mat(spikes(:));
            
            if length(denoised_voltage) ~= length(time)
               str = sprintf('Warning: voltage and time lengths don''t match which may make things a little funky');
               displayErrorMsg(str); 
            end
            
         catch ME
            displayErrorMsg('You''ve fucked something up - check denoiseVoltage:threshold');
         end
                                 
      case 'wavelets'
         try
            % use wavelets to denoise the voltage timeseries
            name   = params.mother_wavelet.value; % mother wavelet
            levels = params.wavelet_level.value;   % level/scale to go up to
            number = params.wavelet_number.value;  % number wavelet to use in family
            method = 'modwtsqtwolog'; % you can add this as a configurable parameter
            [~,filtname] = getMotherWavelet(name, number);
            if strcmpi(method, 'modwtsqtwolog')
               W   = modwt(double(tseries.data), filtname, levels); 
               denoised_voltage = toVec(wden(W, 'modwtsqtwolog', 'h', 'mln', levels, filtname));
            else
               [C, L]  = wavedec(tseries, levels, name); % wavelet decomp - C is time, L is num coeffs by level
               denoised_voltage = wden(C, L, method, 'h', 'mln', levels, name);
            end
         catch ME
            displayErrorMsg('You''ve fucked something up - check denoiseVoltage:wavelets');
         end
         
      otherwise
         str = 'Lazy me - denoising technique not implemented, ignoring';
         displayErrorMsg(str);
         denoised_voltage = [];
   end
   % Print 'done' message
   str = sprintf( '\tDone!\n');
   cprintf( 'Keywords', str );
end