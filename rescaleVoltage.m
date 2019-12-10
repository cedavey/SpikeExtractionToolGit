% v = rescaleVoltage(v, time, dt, params)
%
% Rescale voltage using user specified method from gui

% How flat do you want it?
% variance of peaks changes but still has small peaks
% negative changes differently to positive --> can't have both flat
% rules for removing obvious outliers (e.g. v.large amplitude noise)

% - spread in spike peaks changes but noise var doesn't, whichz means the
%   resistance estimate sits differently with regards to peaks
% - how quickly do you want resistance estimate to change? Fairly stable,
%   or quickly so that all peaks tend towards same amplitude where earlier
%   they didn't

function [rescaled_voltage, Rest, params] = rescaleVoltage(tseries, method, params, varargin)
   % methods: 'Variance', 'Particle filter', 'Recursive least squares', 'Recursive mean'

   % Check if options input variable exist
   if nargin > 3, opts = varargin{1}; else, opts = []; end

   switch lower(method)
      case 'variance'
         rescaled_voltage = rescaleVoltageVariance(tseries, params, opts);
         Rest = [];

      case {'particle filter', 'recursive least squares', 'recursive mean'}
         str = sprintf( '\tRescaling voltage, this may take a whiles...\n' );
         printMessage('off', 'Keywords', str);
         [rescaled_voltage, Rest, ~, params] = rescaleVoltageRecursive(tseries, params, method, opts);

      case 'plot sigma'
         plotVoltageSigma(tseries, params);
         rescaled_voltage = [];
         Rest = [];

      otherwise
         str = 'Nah man, no such rescaling technique';
         displayErrorMsg(str);
         rescaled_voltage = [];
   end
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

function [vrescale, Rest_vec, tpeak_vec, params] = rescaleVoltageRecursive(tseries, params, method, opts)
   % Check that options variable 'opts' contains fields
   if isfield(opts,'debugOption'), debug = opts.debugOption; else, debug = 'semi'; end % debug is now selected from the GUI's context menu (right click). Default is 'semi'.
   if isfield(opts,'loadingWindowOn'), progress_window = opts.loadingWindowOn; else, progress_window = false; end
   if isfield(opts,'rescaleOption'), rescaling_method = opts.rescaleOption; else, rescaling_method = 'at_end';end
   if isfield(opts,'auto_params'), auto_params = opts.auto_params; else, auto_params = false; end

   % extract time details - lambda is a forgetting factor; 1 --> normal
   % mean calculation, 0 --> last sample is the estimate, & somewhere in
   % between gives a recursive update, usually it should be set to ~0.9.
   time   = double( tseries.time ); % valid indexing of time was failing when single
   dt     = double( tseries.dt );
   v      = double( tseries.data );
   nT     = min( length(v), length(time) );
   time   = toVec( time(1:nT) );
   v      = toVec( v(1:nT) );    % for some reason time & v can be 1 sample out
   % Initialize noise and define automatic parameters
   [noisemu, noisesig, initialNoiseSamples]   = initNoiseStdDev(v, debug);      % init noise std to identify spikes
   if auto_params
      % Get auto params
      pp = getAutomaticParams(tseries, noisesig, 'rescale', method, params);
      % Save the select_peaks parameter (positive or negative)
      select_peaks = params.select_peaks.value;
      % If returned value is empty, leave the default parameters untouched
      if ~isempty(pp)
         params = pp;
         params.select_peaks.value = select_peaks; % Restore the value (This is not chosen automatically)
      end
   end

   voltoutlier  = params.voltage_magnitude.value;
   glitchthresh = params.glitch_magnitude.value;
   select_peaks = params.select_peaks.value;
   lambda       = params.forgetting_factor.value;
   jumpahead    = params.jump_ahead.value;
   peakfn       = getPeakFn( select_peaks ); % select just positive peaks, just negative, both or separate
   npoles = 0;     % number of poles in AR model of resistance coeff
   nzeros = 2;     % includes the dc component
   nlags  = max(nzeros, npoles+1);
   jumpahead = round( jumpahead / dt ); % convert jump ahead from time into samples

   tnow = datetime('now');
   str = sprintf("\tStarted rescaling, method %s, spikethresh %i, glitchthresh %i, jump ahead %g s, at %s\n", ...
      method, voltoutlier, glitchthresh, jumpahead*dt, datestr(tnow));
   printMessage('off', 'Keywords', str);

   % Initialise the rest of the values
   [vspike,  tspike]     = initSpikes(time, dt, peakfn(v), noisesig*voltoutlier, peakfn, glitchthresh * noisesig, jumpahead); % identify initial spikes
   [Rest, Rcoeff, Rmu, Rstd, Rcov] = initRegress(tspike-tspike(1), vspike, lambda, npoles, nzeros, debug); % est init resistance
   a = Rcoeff(2:npoles+1); b = [Rcoeff(1) Rcoeff(npoles+2:end)];
   Rinit    = Rest(end); % The initial estimate is what might be causing the initial change in rescalings %TODO
   Rcovprev = Rcov;
   vmu      = mean( vspike );
   vstd     = std( vspike );
   pf       = initpf(Rest, Rstd, npoles); % initialise particle filter with initial resistance
   Rvar     = Rstd^2;

   % Particle filtering equations:
   % - predict:  R(t_sp) = a * R(t_prevsp) * (t_sp - t_prevsp) + e(t)
   % - correct:  P(y(t)) ~ N( R, Rsig )
   % Adaptive regression:
   % - update regression coefficients for R at each timestep

   % tracking noise mean, mse & num samples
   noisethresh = voltoutlier * noisesig;
   noisemse    = (noisesig^2)*initialNoiseSamples;
   noiseN      = initialNoiseSamples; % Used to be 1. Changed to the actual length of the initial white noise (initNoiseStdDev)

   lambda      = ternaryOp( lambda==1, 1-dt, lambda ); % make lambda 1 timestep less than 1
   R_prob      = 0.03;
   remsamp     = nT;
   currt       = 1;
   pprevt       = [1 1];
   vpeak_vec   = [];    % list of all voltage peaks
   tpeak_vec   = [];    % list of time each voltage peak occurred
   Rest_vec    = Rest(end);  % resistance estimate for current timestep
   Rprev       = Rest(end)'; % resistance estimate for prev timestep
   avgpeaks    = 20;    % number of peaks to use when new piece of piecewise regression
   tstartpiece = 0;     % time the current linear piece started (time shift each to 0)
   timeinpiece = 0;     % counts number of samples in linear regression piece
   numpieces   = 1;     % number of pieces in piecewise linear regression
   vrescale = v;

   % go through timeseries, extract voltage peaks & update resistance
   nP = 0;              % num voltage/spike peaks we've processed so far
   % Create new progress window
   if progress_window
      waitbar_handles = waitbar(0,[method ' is 0% done'],...
         'Name','rescaleVoltage - close me to stop');
   end
   prevProg = 1; % Previous progress, start with 1 because it's the first loop. It will be updated every 10 percent
   firstRun = 1;
   secondRun = 0;
   while currt < nT
      % Update Progress window
      if progress_window
         try
            waitbar_handles = waitbar(currt/nT, waitbar_handles,...
               [method ' is ', num2str(round( currt/nT*100 )), '% done']);
         catch
            delete(waitbar_handles);
            str = ['The process has been manually stopped on time: ',num2str(currt*dt),' seconds.\n'];
            printMessage('off','Error',str);
            break;
         end
      end

      % Print current percentage
      if 10*floor(currt/nT*10) ~= prevProg
         str = sprintf( '\t%s is approx %d%% done\n', method, 10*floor(currt/nT*10) );
         printMessage('off','Text', str);
         prevProg = 10*floor(currt/nT*10);
      end

      %% Shifting the 'noise period' until a reasonable number of samples contain white noise
      % find where noise ends - update est of noise mean & std dev
      startsp = find( peakfn(v( currt:end) ) >= noisethresh, 1, 'first' );
      % no more peaks left in timeseries so exit loop
      if isempty( startsp )
         break;
      end
      prevt  = currt;               % copy current time before it's updated
      pprevt = [pprevt(end) prevt];
      currt  = currt + startsp - 1; % update current time to where next spike is
      % if start of spike is 1st index then prevt==currt & we get stuck
      currt  = ternaryOp( prevt==currt, currt+1, currt );

      % update noise mean & variance estimate - keep track of mean
      % because noise has subthreshold spikes so mean may not be 0
      noise = v(prevt:currt-1);
      [noisemu, noisesig, noisemse, noiseN, is_noise, returned_noise] = updateNoiseStats(noise, noisemu, noisesig, noisemse, noiseN, lambda);
      % If the length wasn't at least 30** samples, it will concatenate
      % them to our noise_vector vector until we have 30** samples. It
      % has restored the noise estimate to whatever was before this
      % new noise try. **NOTE: 30 means, min_noise.
      if ~isempty(returned_noise), noise = returned_noise; end
      if ~is_noise
         if ~exist('noise_vector','var')
            noise_vector = noise;
         else
            % If there is already a vector being saved, concatenate the
            % new few samples and check again. If it is still too small,
            % it will again just use the previous noise estimate and
            % exit the loop.
            noise_vector = [noise_vector; noise];
            [noisemu, noisesig, noisemse, noiseN, is_noise] = updateNoiseStats(noise_vector, noisemu, noisesig, noisemse, noiseN, lambda);
            % If it is still smaller than min_noise samples, it will keep
            % the noise estimate and exit the loop. If it happens to be
            % larger than min_noise samples, and is noise, then we will
            % clear the nosie_vector and still exit the loop.
            if is_noise, clear('noise_vector'); end
         end
      end
      %%
      noisethreshprev= noisethresh;  % keep last thresh so we keep spike we stopped at
      noisethresh    = voltoutlier * noisesig;

      % set prev time pt to curr time pt, & advance current time pt to end of spike
      prevt = currt;
      [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh );
      % Flag to re-check spike is within thresholds and it is not too short
      confirm_is_spike = true;
      
      % Use the new glitch thresrhold that depends on spike amplitude,
      % instead of noise.
      if ~exist('newglitchth', 'var') || isempty(newglitchth)
         newglitchth = glitchthresh * noisesig;
         glitch_plot = [newglitchth currt];
      end
      
      % Store newglitchth if we are plotting in debug option
      if ~contains( debug, 'no', 'ignorecase', true )
         glitch_plot = [glitch_plot; [newglitchth currt]];
      else
         clear glitch_plot;
      end
      
      % This has to be checked in a loop cause when one of them happens, the
      % current spike is updated to the next spike, so thresholds and
      % duration have to be checked again.
      while confirm_is_spike
         % Make it false straight away, so it will only reenter here if the
         % spike was too short. Meaning a new spike was chosen and it still
         % needs to check for glitch.
         confirm_is_spike = false;
         % if this peak is crazy big assume it's a glitch so move to next spike
         countLoops = 0; % Checking how often it enters this loop
         while max( peakfn(vsp) ) > newglitchth
            countLoops = countLoops + 1;
            % move on to next large amplitude event
            startsp = find( peakfn(v( currt:end) ) >= noisethresh, 1, 'first' );
            % no more peaks left in timeseries so exit loop
            if isempty( startsp )
               currt = []; % set to empty since no peaks in remainder of signal
               break;
            end
            prevt   = currt + startsp - 1; % update current time to where next spike is
            [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh );

            if auto_params
               if(countLoops > 4)
                  newglitchth = 1.25 * newglitchth;
                  str = sprintf('\tNew glitch Th: %d\n',newglitchth);
                  printMessage('off', 'Text', str);
               end
            else
               if(countLoops > 29)
                  str = sprintf('\tIt looks like you haven''t chosen the best parameters.\n\tThere were at least %d consecutive spikes that were not\n\tunder the glitch threshold.\n',countLoops);
                  printMessage('off', 'Errors', str);
                  break
               end
            end
         end
         % integrity check - if no peaks in current spike jump ahead a bit & try again
         % - this can happen because we've updated noise stats between when we found
         %   the current spike by its amp being greater than noise, and here
         while max( peakfn(vsp) ) < noisethresh
            % find where noise ends
            startsp = find( peakfn(v( currt:end) ) >= noisethresh, 1, 'first' );
            prevt   = currt + startsp - 1; % update current time to where next spike is
            % no more peaks left in timeseries so exit loop
            if isempty( startsp )
               currt = prevt;
               break;
            end
            % get next spike
            [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh );
         end

         % Check for + and - duration of peak to assess wether it should be
         % considered a spike.
         min_pos_time = 7; % Samples. We can make it a user chosen parameter
         min_neg_time = 5; % Samples. We can make it a user chosen parameter
         sp_too_short = isSpikeTooShort(vsp, peakfn, min_pos_time, min_neg_time);
         while sp_too_short
            % move on to next large amplitude event
            startsp = find( peakfn(v( currt:end) ) >= noisethresh, 1, 'first' );
            % no more peaks left in timeseries so exit loop
            if isempty( startsp )
               currt = []; % set to empty since no peaks in remainder of signal
               break;
            end
            prevt   = currt + startsp - 1; % update current time to where next spike is
            [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh );
            sp_too_short = isSpikeTooShort(vsp, peakfn, min_pos_time, min_neg_time);
            confirm_is_spike = true; % So it will re-enter the while loop and check for glitch thresholds again.
         end
      end

      % if we're here & time is invalid then we're at the end of the signal
      if isempty( currt ) || isempty( prevt )
         currt = numel( time ); % advance time to the end (prob not even nec)
         break;
      end

      prevpeaks   = min( avgpeaks, max(1,nP-1) ); % num peaks for avg's & std's
      [vp,pind]   = max( peakfn(vsp) );   % get index of last max peak
      vpeak = vp(end);
      tpeak = prevt + pind(end) - 1;      % get time ind of last max peak
      nP    = nP + 1;

      vpeak_vec(nP,1) = vpeak;
      tpeak_vec(nP,1) = time(tpeak);

      % update voltage std dev
      [vmu, vvar] = recursiveUpdate( vmu, vstd^2, vpeak, lambda );
      vstd = sqrt( vvar );

      % Now that we have a peak, we can change the glitch threshold to
      % depend on the value of the peaks, rather than the noise.
%       newglitchth = 1.2 * max(vpeak_vec);
      newglitchth = 3 * mean(Rest_vec);

      % if we have at least 2 voltage peaks, update state estimate
      if length( tpeak_vec ) > (nlags+1)
         timeinpiece = timeinpiece + 1;

         switch strip( lower( method ) )
            case {'particle filter'}
               % sometimes the particle filter gets stuck if the
               % regression is 0 except for the DC component, so that
               % all particles end up having the same or v.similar value
               if length(unique(pf.Particles))<40
                  pf = initpf(Rest, Rstd, npoles); % initialise particle filter with initial resistance
                  % pf = initpf(Rest, vstd); % initialise particle filter with initial resistance
               end
               % Predict Rest from prev timestep: x(t-1) -> x(t)

               [a, b, Rcov, Rtmp] = recursiveLeastSquares( (tpeak_vec(end-nlags:end) - tstartpiece), ...
                  vpeak_vec(end-nlags:end), ... % Rest, ... %
                  npoles, nzeros, Rcovprev, ...
                  a, b, lambda );
               % Update regression that estimates next resistance
               Rcoeff = [b(1) a(:)' toVec(b(2:end))'];

               predict( pf, Rcoeff, Rest, ( tpeak_vec(nP-nlags:end) - tstartpiece ), npoles, nzeros );

               % Correct Rest from current obs: x(t) -> y(t)
               correct( pf, vpeak, Rstd );
%                correct( pf, Rtmp, Rstd );
               % Use corrected Rest to update regression coeffs
               Rest  = getStateEstimate(pf);
               Rcurr  = Rest(1); % if multiple lagged states get current only
               % If this is the first run, Rest_vec contains the initalized
               % Rest, which seems to be wrong. Fill Rest_vec with Rcurr,
               % which is the first value that seems alright.
               if firstRun == 1, Rest_vec = ones(size(Rest_vec))*Rcurr; firstRun = 0;end

               % if not modelling change in resistance as linear
               % regression then force it to change slowly
%                Rest  = Rprev*lambda + Rcurr*(1-lambda);

               % prepare for next loop
               Rprev = Rest;
               Rcovprev = Rcov;

            case 'recursive least squares'
               % Update regression that predicts mean of spike peak
               [a, b, Rcov, Rest] = recursiveLeastSquares( ( tpeak_vec(end-nlags:end) - tstartpiece ), ...
                  vpeak_vec(end-nlags:end), ...
                  npoles, nzeros, Rcovprev, ...
                  a, b, lambda );
               Rcovprev = Rcov;
               Rcurr    = Rest;
               % Update regression to estimate next resistance
               Rcoeff   = [b(1) a(:)' toVec(b(2:end))'];

            case 'recursive mean'
               Rest     = Rest_vec(end) * lambda  +  vpeak * (1-lambda);
               Rcurr    = Rest;
         end

         % start a new piece of piecewise regression if prob of
         % obtaining voltage peak from resistance is very small, plus
         % the recent run of peak values is not white, indicating
         % a slow drift and the need for a new regression piece (as
         % opposed to a single large resistance value)
         prob_Rest = normpdf( vpeak, Rcurr, vstd );
         if ( prob_Rest < R_prob && ...
               (timeinpiece>=10 && ~iswhite( [Rest_vec(nP-prevpeaks+1:nP-1); Rest_vec(end)] - vpeak_vec(nP-prevpeaks+1:nP), 'bg', 0.15 ) ) ) ...
               || ( Rcurr <= 0 )
            [Rest, Rcoeff, Rmu, Rstd, Rcov] = initRegress(tpeak_vec(nP-prevpeaks+1:nP) - time(currt), ...
               vpeak_vec(nP-prevpeaks+1:nP), lambda, npoles, nzeros, debug);
            Rcurr = Rest(end);
            Rcoeff_vec(nP,:) = Rcoeff;
            % counter to not allow regression reset for a bit else we come back next timestep
            timeinpiece = 0;
            tstartpiece = time(currt);
            numpieces   = numpieces + 1;
            Rvar  = Rstd^2;
            vmu   = mean( vpeak_vec(nP-prevpeaks+1:nP) );
            vstd  =  std( vpeak_vec(nP-prevpeaks+1:nP) );

            if contains( method, 'particle filter' )
               pf = initpf( Rest, Rstd, npoles ); % initialise particle filter with initial resistance
               % pf = initpf( Rest, vstd ); % initialise particle filter with initial resistance
            end
         end
         [Rmu, Rvar] = recursiveUpdate( Rmu, Rvar, Rcurr, lambda );
         Rstd = sqrt(Rvar);
         Rest_vec(nP,1) = Rcurr;

      else % still getting init samples of Rest to enable modelling
         Rest_vec(nP,1) = Rest(end);
      end

      % update vectors
      if contains( debug, 'semi', 'ignorecase', true )
         vstd_vec(nP,1)   = vstd;
         noise_vec(nP,1)  = noisesig;
      end
      if contains( debug, 'full', 'ignorecase', true )
         vstd_vec(nP,1)   = vstd;
         noise_vec(nP,1)  = noisesig;
         Rcov_vec(nP,:,:) = Rcov;
         Rmu_vec(nP,1)    = Rmu;
         Rstd_vec(nP,1)   = Rstd;
         vmu_vec(nP,1)    = vmu;
         if ~strcmpi( strip(lower( method )), 'recursive mean' )
            Rcoeff_vec(nP,:) = Rcoeff;
         end
      end

      % Rescale in semi-real time (every JA). If using it to rescale in
      % semi-real time, we need to save all the Rest_vec and tpeak_vec in
      % order to plot.
      % Rescaling method options: Every JA, between the peaks, every
      % sample, or at the end (original).
      switch rescaling_method
         case 'JA'
            rscltstart = pprevt(1); % Beginning of the period of time to rescale in this loop
            rscltend = pprevt(end); % End of the period of time to rescale in this loop
            if rscltend > length(time), rscltend = length(time); end % Test if we've reached the end of the recording
            if numel(tpeak_vec) > 1
               Rest_vec_realtime = Rest_vec(end -1 : end); % Resistance estimate between the last two found largest spikes
               tpeak_vec_realtime = tpeak_vec(end -1 : end); % Time of the two last found largest spikes
            else
               % If it is the first loop, then there is only 1 estimated
               % resistance
               Rest_vec_realtime = Rest_vec;
               tpeak_vec_realtime = tpeak_vec;
            end
         case 'peaks'
            % Rescale only the voltage between the last two largest peaks.
            % It will rescale from the zero cross before the spike, instead
            % of using the acutal index of the peak, because that would
            % change the shape of the spikes.
            if numel(tpeak_vec) > 1
               rscltstart = find(time == tpeak_vec(end-1),1,'first'); % Beginning of the period of time to rescale in this loop. It is the index of the second last peak
               Rest_vec_realtime = Rest_vec(end - 1 : end); % Resistance estimate between the last two found largest spikes
               tpeak_vec_realtime = tpeak_vec(end - 1 : end); % Time of the two last found largest spikes
            else
               % If it is the first loop, then there is only 1 estimated
               % resistance
               rscltstart = 1;
               Rest_vec_realtime = Rest_vec;
               tpeak_vec_realtime = tpeak_vec;
            end

            % Find the index of the closest positive transition before the
            % second last peak (rscltstart).
            try
               posv = v(pprevt(1): currt)>0; % Get values above 0 (ones), and below (zeros)
            catch
               posv = v(pprevt(1): min(currt, length(v)))>0; % Could ask for this in the previous line and avoid the try-catch statement, but tested with tic toc and it takes longer.
            end
            changePos  = find(diff([0; posv; 0])==1);    % location of change from 0 to 1
            rscltstart = changePos(find(changePos < (rscltstart - pprevt(1)), 1, 'last')); % The actual spike starts on the closest positive transition before the peak
            % Find the index of the closest positive transition after the
            % last peak (rscltend)
            rscltend = find(time == tpeak_vec(end),1,'first'); % End of the period of time to rescale in this loop. It is the index of the second last peak
            rscltend = changePos(find(changePos > (rscltend - pprevt(1)), 1, 'first')); % The actual spike ends on the closest positive transition after the peak
            % We only looked for indexes within the current period, so we
            % have to adjust the indexes by the start of current period.
            rscltstart = rscltstart + pprevt(1);
            rscltend = rscltend + pprevt(1);
            % Check our rescaling period is within the limits of the
            % recording
            if (isempty(rscltstart)||rscltstart < 0), rscltstart = pprevt(1); end
            if (isempty(rscltend) || rscltend > currt), rscltend = currt; end % Test if we've reached the end of the section
      end

      if ~strcmp('at_end',rescaling_method)
         % Rescale
         [vv, ~, ~] = doRescale(v( rscltstart:rscltend ),...
            tpeak_vec_realtime, Rest_vec_realtime,...
            time( rscltstart:rscltend ));
         vrescale(rscltstart:rscltend) = vv;

         % Plot the last rescaled section. (If debug is full)
         if strcmp('full',debug)
            figure(1);
            if numel(tpeak_vec) > 1 && currt < numel(v)
               plot(time(rscltstart:rscltend), vrescale(rscltstart:rscltend));
               hold(gca,'on');
               plot(time(rscltend : currt), v(rscltend : currt),'g');
               plot([time(find(time == tpeak_vec(end-1),1,'first')) time(find(time == tpeak_vec(end),1,'first'))],...
                  [v(find(time == tpeak_vec(end-1),1,'first')) v(find(time == tpeak_vec(end),1,'first'))],'ro','LineWidth',2);
               plot(time(currt), vrescale(currt),'kx','LineWidth',2);
               plot(time(pprevt), vrescale(pprevt),'mx','LineWidth',2);
               plot(time(rscltstart:rscltend), newglitchth*ones(1,1+rscltend-rscltstart));
               xlim([time(max(1,rscltstart - 100)) time(currt)]);
               ylim([-3 3]);
               box off;
               grid off;
               legend({'Rescaled voltage','Original voltage','Last 2 peaks','Current time'})
               hold(gca,'off');
            end
         end
      else
         % If rescaling at the end, we still have to fix the time of the
         % rescaling, so it does not change the shape of a spike. I.e. it
         % should not start rescaling with a different value during a peak,
         % rather during the closest zero cross.
         %
         % Find the index of the closest positive transition before the
         % second last peak (rscltstart).
         try
            posv = v(pprevt(1): currt)>0; % Get values above 0 (ones), and below (zeros)
         catch
            posv = v(pprevt(1): min(currt, length(v)))>0; % Could ask for this in the previous line and avoid the try-catch statement, but tested with tic toc and it takes longer.
         end
         changePos  = find(diff([0; posv; 0])==1);    % location of change from 0 to 1% Find the index of the closest positive transition after the
         % last peak (rscltend)
         rscltend = find(time == tpeak_vec(end),1,'first'); % End of the period of time to rescale in this loop. It is the index of the second last peak
         rscltend = changePos(find(changePos > (rscltend - pprevt(1)), 1, 'first')); % The actual spike ends on the closest positive transition after the peak
         % We only looked for indexes within the current period, so we
         % have to adjust the indexes by the start of current period.
         rscltend = rscltend + pprevt(1);
         % Check our rescaling period is within the limits of the
         % recording
         if (isempty(rscltend) || rscltend > currt), rscltend = currt; end % Test if we've reached the end of the recording
         % Update the index where Rest changes
         tpeak_vec(end) = rscltend * dt;
      end

   end
   % Close progress window
   if progress_window
      try
         close(waitbar_handles);
      catch ME
         if strcmp('Invalid figure handle.',ME.message)
            if exist('waitbar_handles','var')
               clear('waitbar_handles');
            end
         else
            runtimeErrorHandler(ME);
         end
      end
   end

   if ~contains( debug, 'no', 'ignorecase', true )
      figure;
      subplot(221), hold on;
      plot(time,v);
      try
         % If the plot doesn't work at this point is most likely because there
         % were no spikes found.
         if length(vpeak_vec) < 100, marker_size = 8; else, marker_size = 4; end
         plot(tpeak_vec, vpeak_vec, 'kx', 'markersize', marker_size);
         plot(tpeak_vec, Rest_vec,  'm',  'linewidth', 3);
         plot(tpeak_vec, noise_vec * voltoutlier,  'r');
         plot(tpeak_vec, -noise_vec * voltoutlier,  '--r');
         % plot(tpeak_vec, noise_vec * glitchthresh, 'g');
         % plot(tpeak_vec, noise_vec * params.glitch_magnitude.value, 'r');
         plot(glitch_plot(:,2) .* dt, glitch_plot(:,1));
         plot(glitch_plot(:,2) .* dt, -glitch_plot(:,1),'--r');
         % plot(tpeak_vec, -noise_vec * glitchthresh, '--r');
         legend('Voltage', 'V peaks', 'R estimate', 'Th+', 'Th-');
         if numel(tpeak_vec) > 1, xlim( [ tpeak_vec(1) tpeak_vec(end) ] );end
      catch ME
         if strcmp('Vectors must be the same length.', ME.message) && isempty(tpeak_vec)
            close;
            % error('It looks like there were no spikes found. Check the parameter selection.\nThe value of ''currt'' is: %d', currt);
            str = sprintf('It looks like there were no spikes found. Check the parameter selection.\nThe value of ''currt'' is: %d', currt);
            displayErrorMsg(str);
            vrescale = [];
            return
         else
            runtimeErrorHandler(ME);
         end
      end

      subplot(222), hold on;
      plot(tpeak_vec, vpeak_vec, 'k.','markersize',marker_size);
      plot(tpeak_vec, Rest_vec,  'm', 'linewidth', 3);
      plot(tpeak_vec, noise_vec * voltoutlier,  'r');
      % plot(tpeak_vec, noise_vec * glitchthresh, 'r');
      plot(glitch_plot(:,2) .* dt, glitch_plot(:,1), 'r');
      legend('voltage peaks', 'resistance est');
      if numel(tpeak_vec) > 1, xlim( [ tpeak_vec(1) tpeak_vec(end) ] );end

      subplot(223)
      plot(tpeak_vec, noise_vec);
      title('Noise \sigma');
      if numel(tpeak_vec) > 1, xlim( [ tpeak_vec(1) tpeak_vec(end) ] );end

      subplot(224),
      plot(tpeak_vec, vstd_vec);
      title('Voltage \sigma');
      if numel(tpeak_vec) > 1, xlim( [ tpeak_vec(1) tpeak_vec(end) ] );end
   end

   % Moved this bit of code to function 'doRescale' so we can use the
   % function recursively and rescale in semi-real time (every JA).

   % interp doesn't like going outside the input timebounds, so pad Rest
   % with end results & pad time to avoid getting NaNs (changes slowly so ok)
   if strcmp('at_end',rescaling_method)
      [vrescale, Rest_vec, tpeak_vec] = doRescale(v, tpeak_vec, Rest_vec, time);
   else
      %Need to restore Rest to return it properly.
      if tpeak_vec(1) > time(1)
          tpeak_vec = [time(1); tpeak_vec(:)];
          Rest_vec  = [Rinit; Rest_vec(:)];
      end
      if tpeak_vec(end) < time(end)
          tpeak_vec = [tpeak_vec(:); time(end)];
          Rest_vec  = [Rest_vec(:); Rest_vec(end)];
      end
      Rest_vec  = interp1( tpeak_vec, Rest_vec, time );
   end

   tnow = datetime('now');
   str = sprintf("\tFinished rescaling at %s\n", datestr(tnow));
   printMessage('on', 'Keywords', str);

end

% Actually performs the voltage rescaling. The function can be called at
% the end of the resistance estimation or every period, e.g. every jump
% ahead. Returns the rescaled voltage over the specified period.
%
% Inputs:
%  v        - the original voltage or section to be rescaled
% tpek_vec_realtime
%  Rest_vec  - resistance estimate
%  tpeak_vec - time in seconds between which the rescaling is performed. It
%              is expected to be a 2 elements vector. If it is longer it is
%              fine. Rest_vec should be the same size.
%  time      - time vector in seconds. Same size as v.
% Outputs:
%  vrescale - rescaled voltage
%  Rest_vec - Updated Rest vector is an interpolation between the input
%             values of Rest. Can be dismissed.
%  tpek_vec - Updated tpeak_vec is an interpolation between the input
%             values of tpeak_vec. Can be dismissed.
function [vrescale, Rest_vec, tpeak_vec] = doRescale(v, tpeak_vec, Rest_vec, time)
   % interp doesn't like going outside the input timebounds, so pad Rest
   % with end results & pad time to avoid getting NaNs (changes slowly so ok)
   if tpeak_vec(1) > time(1)
      tpeak_vec = [time(1); tpeak_vec(:)];
      Rest_vec  = [Rest_vec(1); Rest_vec(:)];
   end
   if tpeak_vec(end) < time(end)
      tpeak_vec = [tpeak_vec(:); time(end)];
      Rest_vec  = [Rest_vec(:); Rest_vec(end)];
   end

   % Remove zeros
   zidx = tpeak_vec == 0;
   tpeak_vec(zidx) = [];
   Rest_vec(zidx) = [];

   Rest_vec  = interp1( tpeak_vec, Rest_vec, time );
   vrescale  = v(:) ./ Rest_vec(:);
end

% Get the current spike. Or current bundle of samples that begins with a
% sample greater than the noise threshold. Get at least 'jumpahead'
% samples, but end on a noise sample. So if there's a run of large
% amplitude samples don't stop half way through, but get 'em all.
% Note: the function assumes that we're at the start of the spike, or large
%       amplitude event
% Inputs:
%  v           - voltage trace
%  currt       - index of where we're up to in the voltage trace
%  jumpahead   - size of spike to get
%  noisethresh - below this threshold is noise, above is something else
% Outputs:
%  vsp         - vector of 'spike'
%  currt       - get next current time, set to sample just after last spike
%                sample
function [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh )
   % we're at next spike - only want to model peaks, so just
   % get next 20 or so samples, which should contain the whole spike,
   % and then from there find the next noise sample - this way we get a
   % small set of samples that should contain the whole spike, and then
   % the next time through the loop we'll be starting in noise, and we
   % can compare the set of samples against current resistance & do PF etc
   % (it also makes sure wacky bits of noise only impact us once!)

   remsamp = length(v) - prevt;  % remaining num samples
   jump    = double(min( remsamp, jumpahead )); % When jump was 'single' instead of 'double' it was causing an error if it was used as an index

   % prevt   = prevt + find( peakfn(v(prevt+jump:end)) < noisethresh, 1, 'first');
   tnoise  = find( peakfn(v(prevt+jump:end)) < noisethresh, 1, 'first');

   % if there's no noise before the end of the sample, jump to the end
   jump    = jump + ternaryOp( isempty( tnoise ), numel(v(prevt+jump:end)), tnoise );
   endsp   = jump - 1;
   vsp     = v(prevt : prevt+endsp);    % exract current spike
   currt   = prevt + endsp + 1;
end

% Checks if the duration of positive and negative parts of the highest
% spike are reasonable.
%
% Inputs:
%   vsp       - vector of 'spike'
%   peakfn    - peak function
% Outputs:
%   too_short - flag is true if spike is too short on either sign (+ or -)
function too_short = isSpikeTooShort(vsp, peakfn, min_pos_time, min_neg_time)
   pos_time = min_pos_time;
   neg_time = min_neg_time;

   [~, idx] = max( peakfn(vsp) );  % Get the highest peak
   posv = vsp>0; % Get values above 0 (ones), and below (zeros)
   % find where voltage goes from -ve to +ve
   changePos  = find(diff([0; posv; 0])==1);    % location of change from 0 to 1
   stidx = changePos(find(changePos < idx, 1, 'last')); % The actual spike starts on the closest positive transition before the peak
   if (isempty(stidx)||stidx < 0), stidx = 1; end
   endidx = changePos(find(changePos > idx, 1, 'first')); % The actual spike ends on the closest positive transition after the peak
   if (isempty(endidx)||endidx > length(vsp)), endidx = length(vsp); end
   % too_short is true if there is less consecutive 1's than the desired
   % length of the peak (positive part of the spike), or less consecutive
   % 0's than the desired length of the valley (negative part of the
   % spike)
   too_short = (sum(posv(stidx:endidx)) < pos_time) || (sum(~posv(stidx:endidx)) < neg_time);
end

%--- loop through noise & update mean & mse recursively ---
function [mu_curr, sig_curr, mse_curr, N_curr, is_noise, returned_noise] = updateNoiseStats(new_noise, mu_prev, sig_prev, mse_prev, N_prev, lambda)
   % Check if "new_noise" is actually white noise, adjust the length of the
   % samples like we did in initNoiseStdDev (actually using same function).
   % If can't find at least 30 samples (or the length of the new_noise
   % signal), then use the whole new_noise. See what to do in that case.

   % If new_noise contains too few samples, mantain the previous noise
   % value.
   is_noise = true;
   min_noise = 30;
   returned_noise = [];
   if length(new_noise) < min_noise
%       str = sprintf('\tCouldn''t find %d samples between Spikes\n',min_noise);
%       printMessage('off','Text',str);
      is_noise = false;
      mu_curr = mu_prev;
      sig_curr = sig_prev;
      mse_curr = mse_prev;
      N_curr = N_prev;
      return
   else
      try
         [~, ~, ~, nn] = initNoiseStdDev(new_noise, 'none', min_noise, 0.15);
         new_noise = nn; % If it didn't trigger the catch statement
      catch E
         if contains(E.message, 'find a single white period')
            % It means it didn't find 30 samples (or whatever number) of
            % white noise. If it came here, we will not care and just treat
            % the whole signal new_noise as if it was actual noise.
%             str = sprintf('\tCouldn''t find at least %d samples of white noise. %d \n', min_noise, N_prev + length(new_noise));
%             printMessage('off','Errors',str);

            is_noise = false;
            mu_curr = mu_prev;
            sig_curr = sig_prev;
            mse_curr = mse_prev;
            N_curr = N_prev;
            return
%             is_noise = false;
%             mu_curr = mu_prev;
%             sig_curr = sig_prev;
%             mse_curr = mse_prev;
%             N_curr = N_prev;
%             idx = randi( numel(new_noise) - floor(min_noise/3),1 );
%             returned_noise = new_noise(idx: idx + floor(min_noise/3) - 1);
%             return
            % Try again with a smaller no of samples
%             [~, ~, ~, nn] = initNoiseStdDev(new_noise, 'none', 10, 0.15);
         else
            runtimeErrorHandler(E);
         end
      end
   end

   recursive = true;
   N_curr    = N_prev + length( new_noise );

   % update noise using either recursive mean or, if assumed stationary,
   % using normal mean/sig calculations in real time
   if recursive
      % ASSUMES NOISE IS NON-STATIONARY
      % recursive updating of noise, with forgetting
      % init curr vars cuz I seem to have a problem with new noise
      % sometimes being empty?!
      mu_curr = mu_prev; sig_curr = sig_prev; mse_curr = mse_prev;
      for ni=1:length(new_noise)
         % recursive noise estimates
         [mu_curr, var_curr] = recursiveUpdate( mu_prev, sig_prev^2, new_noise(ni), lambda );
         mse_curr = mse_prev + var_curr;
         mse_curr = mse_prev + (new_noise(ni) - mu_curr) * (new_noise(ni) - mu_prev);
         mu_prev  = mu_curr;
         sig_prev = sqrt(var_curr);
         mse_prev = mse_curr;
      end
   else
      % ASSUMES NOISE IS STATIONARY
      % normal mean/sig estimates in real time
      mu_curr  = ( mu_prev * N_prev + sum( new_noise ) ) / N_curr;
      mse_curr = mse_prev + (sum(new_noise) - mu_curr)*(sum(new_noise) - mu_prev); % mse_curr = mse_prev + sum( (new_noise - mu_curr) .* (new_noise - mu_prev) );
   end

   sig_curr = sqrt( mse_curr / (N_curr - 1) );
end

% Initialise particle filter object. If a pf is input then reset it
function pf = initpf(Rest, Rsig, npoles, pf)
   % particle matrix has 1st column the most recent timestep, and later
   % columns as earlier timesteps, so reverse Rest to comply
   % using PF to evolve to next step so current step becomes last timestep,
   % and then we don't need the very oldest cuz it falls of the end
   Rest = Rest(end:-1:1);
   Rest = Rest(1:npoles+1); % +1 for current timestep (& previous npoles)
   if ~exist('pf','var') || isempty(pf)
      % min number of particles recommended is 1000, unless processing time is impt
      numParticles = 1e3;
      pf = robotics.ParticleFilter; % initiate a new particle filter
   end
   npoles = length( Rest );
   % mean for each lagged timestep is given by Rest, i.e. particles ~ N( Rest, Rsig )
   initialize(pf, numParticles, ternaryOp( Rest>0, Rest, 0), ones(size(Rest))*Rsig, 'StateOrientation', 'row');

   if nargin<4
      % object estimates state from 'mean' of particles, or particle with 'maxweight'
      pf.StateEstimationMethod    = 'mean'; % if maxweight no cov returned from pf predict
      % resampling method options: 'multinomial', 'residual', 'stratified', 'systematic'
%       pf.ResamplingMethod         = 'systematic'; % no idea what this is?!
      %  evolve particles to get new estimates of the resistance
      pf.StateTransitionFcn       = @evolveResistance;
      % obs likelihood - likelihood of each particle state given new obs
      pf.MeasurementLikelihoodFcn = @noiseLikelihood;
   end
end

%  initialise with regression on first few samples - else we might miss
%  the first few spikes while the recursive algorithm is converging
function [Rest, Rcoeff, Rmu, Rsig, Rcov] = initRegress(tspike, vspike, lambda, npoles, nzeros, debug, Rcov)
   dorecursive = true; % can assume stationary for init section, or do recursively

   % hopefully this will never happen, but just in case!
   if length(vspike)==1
      Rcoeff(2) = 0;
      Rcoeff(1) = vspike;
      Rsig = vspike;
      Rmu  = vspike;
      Rcov = eye(2,2) / Rsig;
      Rest = vspike;
      str  = sprintf( '\tInit regression for R estimate only has 1 sample!\n' );
      printMessage('off', 'Errors', str );
      return;
   end

   % get init estimate of VAR model to initialise recursive least squares
   p      = max(npoles,nzeros);      ndim = npoles + nzeros;
   if numel(tspike) > 4*p+2
      % Validate length of the time series (even if it was larger than 1)
      [Rcoeff, ~, vest] = getARcoeff_2vars( vspike, tspike, p );
   else
      Rcoeff(2) = 0;
      Rcoeff(1) = vspike(end);
      Rsig = vspike(end);
      Rmu  = vspike(end);
      Rcov = eye(2,2) / Rsig;
      Rest = vspike(end);
      str  = sprintf( '\tInit regression for R estimate has too few samples!\n' );
      printMessage('off', 'Errors', str );
      return;
   end

   a      =  Rcoeff(2:npoles+1)'; % autoregression coeffs
   b      = [Rcoeff(1) Rcoeff(p+2:(p+nzeros))']; % DC + input coeffs

   % if the estimate is shit (e.g. negative) then redo as univar AR &
   % regression (i.e. separate out from VAR into separate equations)
   X = makeRecursiveLeastSquaresMatrix( tspike, vspike, length(vspike), npoles, nzeros );
   Rest = [b(1) a, b(2:end)'] * X;
   if Rest <= 0 || Rest < mean(vspike) - 2*std(vspike) ||  Rest > mean(vspike) + 2*std(vspike)
      a = getARcoeff( vspike, npoles );
      b = a(1); a = a(2:end)'; % get rid of dc

      % regress time onto spikes
      if nzeros>1
         clear X;
         for ti=1:(nzeros-1)
            X(:,ti) = tspike((nzeros-ti+1):end-ti+1);
         end
         c = regress( vspike(nzeros:end), X );
         b = [b c];
      end
   end
   if nargin<7 || isempty(Rcov)
      Rcov = eye( ndim ) / std(vspike); % inverse cov of predictor matrix
      Rcov = eye( ndim ); % inverse cov of predictor matrix
   end
   nT     = length(vspike);

   Rmu    = mean(vspike); Rmu_prev  = Rmu;
   Rvar   = var(vspike);  Rvar_prev = Rvar;
   Rest   = vest(end-npoles:end); % init to something just in case insuffic samples
   Rsig   = std(vspike);

   if ~strcmp('none',debug)
      Rmu_vec   = zeros(nT,1);       Rvar_vec  = zeros(nT,1);
      a_vec     = zeros(nT,npoles);  b_vec     = zeros(nT,nzeros);
      Rest_vec  = zeros(nT,1);       Rcov_vec  = zeros(nT,ndim,ndim);
      a_vec(p,:)= a;                 b_vec(p,:)= b;
   end
   if dorecursive
      for ti=(p+1):nT
         % recursive spike mean/std estimates
         [Rmu, Rvar] = recursiveUpdate( Rmu, Rvar, vspike(ti), lambda );

         % linear regression for spike std
         [a, b, Rcov, Rest] = recursiveLeastSquares( ( tspike(ti-p:ti) - tspike(1) ), ...
            vspike(ti-p:ti), ...
            npoles, nzeros, Rcov, ...
            a, b, 1 );

         if ~strcmp('none',debug)
            Rmu_vec(ti,1)    = Rmu;
            Rvar_vec(ti,1)   = Rvar;
            a_vec(ti,:)      = a(:)';
            b_vec(ti,:)      = b(:)';
            Rcov_vec(ti,:,:) = Rcov;
         end
         Rest_vec(ti,:) = Rest;
      end
      Rest   = Rest_vec(end-p:end); % return num lagged inputs
      % Rest = Rest_vec(p+1:2*p+1); % want beginning of valid R estimates, not 400 seconds in or whatever
      Rcoeff = [b(1) a(:)' toVec(b(2:end))'];
      Rsig   = ternaryOp( Rvar<0, std(vspike), sqrt( Rvar ) );
   end
end

% get initial spikes to kickstart regression with - gets the smallest subset
% of voltage samples that is required to extract a specified number of spikes
function [vspike, tspike] = initSpikes(time, dt, v, noiseoutlier, peakfn, glitchthresh, jumpahead)
   % num samples to start with - double until we have enough spikes
   % - if we assume 30 spikes/s to start with, then calc num samples
   %   (we expand time if the actual spiking rate isn't high enough, and
   %    some of these peaks will be from the same spike, so make rate kinda high)
   initp = 1e2;  % number of peaks we require (will drop when remove multiple peaks per spike)

   rate  = 30;
   nS    = round( initp / rate / dt );
   maxS  = round( length(time) * 0.4 ); % don't wanna use too much of v to get spikes!
   fail  = true;
   % for noiseoutlier threshold, find 'initp' voltage peaks that are
   % greater than the threshold, & use these to initialise spikes
   % - but if we can't find 'initp' peaks within a certain proportion of
   % the voltage timeseries, reduce the threshold
   while fail
      while (sum( peakfn(v(1:nS)) > noiseoutlier & peakfn(v(1:nS)) < glitchthresh ) < initp) && nS<maxS
         nS   = min(maxS, nS*2 );
      end
      % if num samples (i.e. nS) is smaller than max then we've reached the
      % target number of initial peaks
      if nS < (maxS+1)
         fail = false;
         % if we exited while loop coz we ran out of timeseries, reduce thresh
      else
         noiseoutlier = noiseoutlier * 0.9;
      end
   end
   vtmp    = v(1:nS); ttmp = time(1:nS);
   isglitch= peakfn( vtmp ) > glitchthresh;
   vtmp    = vtmp( ~isglitch );
   ttmp    = ttmp( ~isglitch );
   isspike = peakfn(vtmp) > noiseoutlier;
   vtmp    = vtmp(isspike);
   ttmp    = ttmp(isspike);

   % extract just the spike peaks - for each run of voltages that are
   % pretty close in time, get just the largest in the run so that we're
   % only getting one sample per spike
   sptime = dt*200; % peaks within this period are assumed part of same spike
   vspike = []; tspike = [];
   while ~isempty(vtmp)
      tdiff = diff(ttmp);
      %
      endsp = find( tdiff > sptime, 1, 'first' );
      % get voltage for this spike only - if no timesteps > dt -> 1 spike left
      if isempty( endsp ), endsp = length( vtmp ); end
      vsp = peakfn( vtmp(1:endsp) );

      % get peak voltage of this spike & time of it's occurrence
      [~, maxi]       = max( vsp );
      vspike(end+1,1) = vsp(maxi);
      tspike(end+1,1) = ttmp(maxi); % index ok for time too coz vsp started at vtmp(1)

      %       % get peak negative voltage
      %       [~, maxi] = max( -vsp( vsp < 0 ) );
      %       if ~isempty(maxi)
      %          vspike(end+1,1) = vsp(maxi);
      %          tspike(end+1,1) = ttmp(maxi); % index ok for time too coz vsp started at vtmp(1)
      %       end

      % ditch spike we just processed
      vtmp(1:endsp) = [];
      ttmp(1:endsp) = [];
   end

   %    % jumpahead param tells us how far ahead we look to get the next peak
   %    % from - e.g. if it's 10 seconds, then we get max voltage every 10
   %    % seconds --> find out how many seconds we've got here & get the max
   %    % voltages from it (jumpahead is in samples, not time)
   %    vspike = zeros(initp, 1); tspike = zeros(initp,1);
   %    starttime = 1;
   %    jumpahead = floor( jumpahead / 4 ); % don't wanna end too far into to recording
   %    for ns = 1:initp
   %       endtime   = starttime + jumpahead - 1;
   %       vsegment  = v(starttime:endtime);
   %       % get rid of voltage that's crazy large - i.e. it's a glitch
   %       vsegment( peakfn(vsegment) > glitchthresh ) = 0;
   %       [vspike(ns), i] = max( vsegment );
   %       tspike(ns)  = time( starttime + i );
   %       starttime = endtime + 1;
   %    end
end


% initialise our estimate of noise std dev - assumes amplitudes change
% slowly so that initial samples will be stationary, so just search for the
% longest subset of samples that are white, & estimate noise from that
%
% Inputs:
%  v     - input voltage vector (whatever number of samples of a recording)
%  debug - (optional) if true, prints the number of samples considered for
%           the white noise. It also plots the white noise section.
%  minimumSamples - (optional)
%  alpha - (optional)
%
% Outputs:
%  noisemu - Mean of the noise
%  noisesig - standard deviation of the noise
%  nS - Number of samples found to be white noise
%  noiseVector - Returns the actual vector found to be white noise
%
% Last modified - Artemio 05-Jun-2019
function [noisemu, noisesig, nS, noiseVector] = initNoiseStdDev(v, varargin)
   % get initial noise estimate by finding a white time series at the
   % beginning of the voltage timeseries
   if numel(varargin) > 0, debug = varargin{1}; else, debug = 'semi'; end
   if numel(varargin) > 1, minimumSamples = varargin{2}; else, minimumSamples = 100; end
   if numel(varargin) > 2, alpha = varargin{3}; else, alpha = 0.05; end

   nS = min(minimumSamples, (length(v) - 1)); % Initial num of samples
   firstS = 1; % Subset starts at this sample
   vtmp = v(firstS : firstS + nS);

   if numel(vtmp) == 1
      % input only has 1 sample. Can't check if it's white noise.
      noisemu = vtmp;
      noisesig = 0;
      nS = 1;
      noiseVector = vtmp;
      return
   end

   while ~iswhite(vtmp, 'lb', alpha, 3)
      % Start with a short number of samples (i.e. 30), if it doesn't
      % evaluate as not-white, shift the window 30 samples ahead.
      firstS = firstS + floor(nS/3);
      try
         vtmp = v(firstS : firstS + nS);
      catch E
         if contains(lower(E.message), 'index exceeds the number of array elements')
            error(['Couldn''t find a single white period of ' num2str(nS) ' continuous samples in the recording.'])
         else
            runtimeErrorHandler(E,'rethrow');
         end
      end
   end

   % If it is white noise, double the number of samples to get the highest
   % possible number of samples that are white noise.
   nSdouble = nS;
   while iswhite(vtmp)
      % If it reached this point, it means that doubling the samples in the
      % previous loop worked and the timeseries kept being white noise,
      % hence it is safe to update nS.
      nS = nSdouble;
      nSdouble = nS * 2; % Doubling the samples to test if still white
      if firstS + nSdouble > length(v)
         % If doubling the samples exceeds the size of the time series,
         % stop here. nS hasn't been updated, so it's fine.
         break;
      end
      vtmp = v(firstS : firstS + nSdouble);
      % If vtmp with double samples is still white noise, the while loop
      % will continue and we'll update nS. If the new vtmp is not white
      % noise, then the while loop will stop here. We haven't updated nS.
   end

   if strcmp('full', debug)
      figure('Name', 'Time series white noise initial period');
      plot(v(firstS : firstS + nS));
      drawnow;
      ylim([-0.15 0.15]);
      str = sprintf('\tWhite noise initialization included %d samples, starting in sample number %d \n', nS, firstS);
      printMessage('off','Text', str);
   end

   noisesig =  std( v(firstS : firstS + nS) );
   noisemu  = mean( v(firstS : firstS + nS) );
   noiseVector = v(firstS : firstS + nS);
end

function rescaled_voltage = rescaleVoltageVariance(tseries, params, opts)
   rescaled_voltage = [];
   time = tseries.time;
   dt   = tseries.dt;
   avg_window  = params.moving_avg_window.value; % size of avg window to incl in calc (seconds)
   skip_window = params.skip_window.value;       % how much to skip ahead for next calc (seconds)
   posneg      = params.separate_pos_neg.value;  % consider +/- voltage separately
   outlier     = params.outlier.value;           % remove volt larger than this is var calc
   dist        = params.distribution.value;      % dist to rescale with - gauss or stud-t
   switch dist
      case {'gauss', 'gaussian', 'norm', 'normal'}
         dist = @std;
      case {'stud-t', 'student''s-t'}
         dist = @studt;
   end
   window      = round(avg_window/dt);     % number of samples to incl in avg
   window      = ternaryOp(iseven(window), window, window+1); % make window even
   skip        = round(skip_window/dt);    % number of samples to incl in skip
   skip        = ternaryOp( iseven(skip), skip, skip+1); % easier if even

   % remove outlier voltages, but retain the values so we can put
   % them back in again later
   voltage     = tseries.data;
   vstd        = std(voltage);
   vthresh     = outlier * vstd;
   v_outInd    = abs(voltage) > vthresh;     % get index of outlier
   %          v_outVal    = voltage(v_outInd);          % record outlier's value
   v_outSign   = sign( voltage(v_outInd) );  % record if outlier's +ve or -ve
   % to avoid 0s when averaging, replace outliers with max value
   voltage( v_outInd ) = vthresh * v_outSign;

   % if posneg is true then smooth +ve & -ve signal components separately
   if posneg
      % try smoothing pos & neg components separately
      posind     = voltage >= 0;
      negind     = voltage <  0;
      vpos       = voltage(posind);
      % make a time vector for positive & negative components
      timepos    = time(posind);
      tpos_      = timepos(ceil(window/2):skip:end-floor(window/2));
      vneg       = voltage(negind);
      timeneg    = time(negind);
      tneg_      = timeneg(ceil(window/2):skip:end-floor(window/2));
      if window ~= 0
         sigpos_ = movingAverageVariance(vpos, window, skip, dist);
         signeg_ = movingAverageVariance(vneg, window, skip, dist);
      else
         sigpos_ = std(vpos);
         signeg_ = std(vneg);
         tpos_   = timepos;
         tneg_   = timeneg;
      end
      % interpolate sigma so it's got same number of samples as voltage
      sigpos   = interp1(tpos_, sigpos_, timepos);
      % interp1 does not like having to extrapolate if times of smaller &
      % larger vectors do not match at beginning and end
      if any(isnan(sigpos))
         tpos_(1) = timepos(1); tpos_(end) = timepos(end);
         sigpos   = interp1(tpos_, sigpos_, timepos);
      end
      signeg   = interp1(tneg_, signeg_, timeneg);
      if any(isnan(signeg))
         tneg_(1) = timeneg(1); tneg_(end) = timeneg(end);
         signeg   = interp1(tneg_, signeg_, timeneg);
      end
      vpos     = tseries.data(posind) ./ sigpos;
      vneg     = tseries.data(negind) ./ signeg;
      rescaled_voltage(posind,1)= vpos;
      rescaled_voltage(negind,1)= vneg;

   else
      % calculate mvg avg std dev at required lowpass frequency
      if window ~= 0
         sigma   = movingAverageVariance(voltage, window, skip, dist);
         % get time vector for downsampled timeseries
         tsmooth = time(ceil(window/2):skip:end-floor(window/2));
      else
         sigma   = nanstd(voltage);
         tsmooth = time;
      end
      % interpolate sigma so it's got same number of samples as voltage
      sig_ma   = interp1( tsmooth, sigma, time );
      % interp1 does not like having to extrapolate if times of smaller &
      % larger vectors do not match at beginning and end
      if sum(isnan(sig_ma))>0
         sig_ma = interp1( tsmooth/(tsmooth(end)-tsmooth(1))*(time(end)-time(1)), sigma, time );
         sig_ma(isnan(sig_ma)) = nanmean(sigma);
      end
      % rescale voltage using sigma
      rescaled_voltage = tseries.data ./ sig_ma;
   end
end

%%
function plotVoltageSigma(tseries, params)

   % voltage_prob = params.voltage_prob.value;
   % inclneg      = params.include_neg_peaks.value;
   % lambda       = params.forgetting_factor.value;

   voltoutlier  = params.voltage_magnitude.value;
   glitchthresh = params.glitch_magnitude.value;
   select_peaks = params.select_peaks.value;
   lambda       = params.forgetting_factor.value;
   % select just positive peaks, just negative, or both
   peakfn       = getPeakFn( select_peaks );
   % % get prob of voltage assuming unit std dev, & scale by noise std dev later
   % if inclneg
   %   voltoutlier  = norminv( 1 - voltage_prob/2, 0, 1 ); % will multiply by std dev later
   % else
   %   voltoutlier  = norminv( 1 - voltage_prob, 0, 1 ); % will multiply by std dev later
   % end
   % extract time details - lambda is a forgetting factor; 1 --> normal
   % mean calculation, 0 --> last sample is the estimate, & somewhere in
   % between gives a recursive update, usually it should be set to ~0.9.
   time   = tseries.time;
   dt     = tseries.dt;
   v      = double( tseries.data );
   nT     = length(v);
   npoles = 0;     % number of poles in AR model of resistance coeff
   nzeros = 2;     % includes the dc component
   nlags  = max(nzeros, npoles+1);

   % Initialise
   [noisemu, noisesig]   = initNoiseStdDev(v);      % init noise std to identify spikes
   [vspike,  tspike]     = initSpikes(time, dt, peakfn(v), noisesig*voltoutlier, peakfn, glitchthresh * noisesig); % identify initial spikes


   debug       = 'semi'; % 'none', 'semi', 'full'

   % tracking noise mean, mse & num samples
   noisethresh = voltoutlier * noisesig;
   noisemse    = sqrt(noisesig);
   noiseN      = 1;

   lambda      = ternaryOp( lambda==1, 1-dt, lambda ); % make lambda 1 timestep less than 1
   R_prob      = 0.01;
   jumpahead   = 5e4; % 5e2; % 5e3; % 24;
   remsamp     = nT;
   currt       = 1;
   vpeak_vec   = [];    % list of all voltage peaks
   tpeak_vec   = [];    % list of time each voltage peak occurred
   avgpeaks    = 20;    % number of peaks to use when new piece of piecewise regression
   timeinpiece = 0;     % number of samples in linear regression piece
   tstartpiece = 0;     % time the current linear piece started
   numpieces   = 1;     % number of pieces in piecewise linear regression


   tnow = datetime('now');
   str = sprintf("Started plotting sigma, spikethresh %i, glitchthresh %i at %s\n", voltoutlier, glitchthresh, datestr(tnow));
   printMessage('off', 'Keywords', str)


   % go through timeseries, extract voltage peaks & update resistance
   nP = 0;              % num voltage/spike peaks we've processed so far
   while currt < nT
      if mod( nP, 1e3 ) == 0
         str = sprintf( '\tPS is approx %d%% done\n', round( currt/nT*100 ) );
         printMessage('off','Text', str);
      end
      % find where noise ends - update est of noise mean & std dev
      startsp = find( peakfn(v( currt:end) ) >= noisethresh, 1, 'first' );
      % no more peaks left in timeseries so exit loop
      if isempty( startsp )
         break;
      end
      prevt  = currt;               % copy current time before it's updated
      currt  = currt + startsp - 1; % update current time to where next spike is
      % if start of spike is 1st index then prevt==currt & we get stuck
      currt  = ternaryOp( prevt==currt, currt+1, currt );

      % update noise mean & variance estimate - keep track of mean
      % because noise has subthreshold spikes so mean may not be 0
      noise = v(prevt:currt-1);
      [noisemu, noisesig, noisemse, noiseN] = updateNoiseStats(noise, noisemu, noisesig, noisemse, noiseN, lambda);
      noisethreshprev= noisethresh;  % keep last thresh so we keep spike we stopped at
      noisethresh    = voltoutlier * noisesig;

      % set prev time pt to curr time pt, & advance current time pt to end of spike
      prevt = currt;
      [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh );
      % if this peak is crazy big assume it's a glitch so move to next spike
      while max( peakfn(vsp) ) > glitchthresh * noisesig
         % move on to next large amplitude event
         startsp = find( peakfn(v( currt:end) ) >= noisethresh, 1, 'first' );
         % no more peaks left in timeseries so exit loop
         if isempty( startsp )
            currt = prevt;
            break;
         end
         prevt   = currt + startsp - 1; % update current time to where next spike is
         [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh );
      end

      % integrity check - if no peaks in current spike jump ahead a bit & try again
      % - this can happen because we've updated noise stats between when we found
      %   the current spike by its amp being greater than noise, and here
      while max( peakfn(vsp) ) < noisethresh
         % find where noise ends
         startsp = find( peakfn(v( currt:end) ) >= noisethresh, 1, 'first' );
         prevt   = currt + startsp - 1; % update current time to where next spike is
         % no more peaks left in timeseries so exit loop
         if isempty( startsp )
            currt = prevt;
            break;
         end
         % get next spike
         [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh );
      end

      prevpeaks   = min( avgpeaks, max(1,nP-1) ); % num peaks for avg's & std's
      [vpeak,pind]= max( peakfn(vsp) );   % get index of peak
      tpeak = prevt + pind - 1;           % get time ind of peak
      nP    = nP + 1;

      vpeak_vec(nP,1) = vpeak;
      tpeak_vec(nP,1) = time(tpeak);

      % if we have at least 2 voltage peaks, update state estimate
      if length( tpeak_vec ) > (nlags+1)
         timeinpiece = timeinpiece + 1;
      end


      % start a new piece of piecewise regression if prob of
      % obtaining voltage peak from resistance is very small, plus
      % the recent run of peak values is not white, indicating
      % a slow drift and the need for a new regression piece (as
      % opposed to a single large resistance value)

      % update vectors

      if contains( debug, 'semi', 'ignorecase', true )
         noise_vec(nP,1)  = noisesig;
      end

   end

   if ~contains( debug, 'no', 'ignorecase', true )
      figure;
      subplot(221), hold on;
      plot(time,v);
      plot(tpeak_vec, noise_vec * glitchthresh, 'r');
      plot(tpeak_vec, noise_vec * voltoutlier, 'g');
      subplot(223)
      plot(tpeak_vec, noise_vec);
      title('Noise \sigma');

   end

   tnow = datetime('now');
   str = sprintf("Finished plotting sigma at %s\n", datestr(tnow));
   printMessage('off', 'Keywords', str);
end
