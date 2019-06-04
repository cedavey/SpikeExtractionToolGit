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

%% TO DO:
% - recursive update of noise in updateNoiseStats can be sped up by
% calcuting new mean & std & weighting by num samples
function [rescaled_voltage, Rest] = rescaleVoltage(tseries, method, params)
% methods: 'Variance', 'Particle filter', 'Recursive least squares', 'Recursive mean'

switch lower(method)
    case 'variance'
        rescaled_voltage = rescaleVoltageVariance(tseries, params);
        Rest = [];
        
    case {'particle filter', 'recursive least squares', 'recursive mean'}
        str = sprintf( '\tRescaling voltage, this may take a whiles...\n' );
        cprintf( 'Keywords', str );
        [rescaled_voltage, Rest] = rescaleVoltageRecursive(tseries, params, method);
        
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

function [vrescale, Rest_vec, tpeak_vec] = rescaleVoltageRecursive(tseries, params, method)
voltoutlier  = params.voltage_magnitude.value;
glitchthresh = params.glitch_magnitude.value;
select_peaks = params.select_peaks.value;
lambda       = params.forgetting_factor.value;
jumpahead    = params.jump_ahead.value;
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
time   = double( tseries.time ); % valid indexing of time was failing when single
dt     = double( tseries.dt );
v      = double( tseries.data );
nT     = min( length(v), length(time) );
v      = toVec( v(1:nT) );    % for some reason time & v can be 1 sample out
time   = toVec( time(1:nT) );
npoles = 0;     % number of poles in AR model of resistance coeff
nzeros = 2;     % includes the dc component
nlags  = max(nzeros, npoles+1);

tnow = datetime('now');
str = sprintf("\tStarted rescaling, method %s, spikethresh %i, glitchthresh %i, jump ahead %g, at %s\n", ...
    method, voltoutlier, glitchthresh, jumpahead, datestr(tnow));
cprintf( 'Keywords', str);


jumpahead = round( jumpahead / dt ); % convert jump ahead from time into samples

% Initialise
[noisemu, noisesig]   = initNoiseStdDev(v);      % init noise std to identify spikes
[vspike,  tspike]     = initSpikes(time, dt, peakfn(v), noisesig*voltoutlier, peakfn, glitchthresh * noisesig, jumpahead); % identify initial spikes
[Rest, Rcoeff, Rmu, Rstd, Rcov] = initRegress(tspike-tspike(1), vspike, lambda, npoles, nzeros); % est init resistance
a = Rcoeff(2:npoles+1); b = [Rcoeff(1) Rcoeff(npoles+2:end)];
Rinit    = Rest(1); % Rest(end);
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

debug       = 'semi'; % 'none', 'semi', 'full'

% tracking noise mean, mse & num samples
noisethresh = voltoutlier * noisesig;
noisemse    = noisesig^2;
noiseN      = 1;

lambda      = ternaryOp( lambda==1, 1-dt, lambda ); % make lambda 1 timestep less than 1
R_prob      = 0.01;
remsamp     = nT;
currt       = 1;
vpeak_vec   = [];    % list of all voltage peaks
tpeak_vec   = [];    % list of time each voltage peak occurred
Rest_vec    = Rest(end);  % resistance estimate for current timestep
Rprev       = Rest(end)'; % resistance estimate for prev timestep
avgpeaks    = 20;    % number of peaks to use when new piece of piecewise regression
tstartpiece = 0;     % time the current linear piece started (time shift each to 0)
timeinpiece = 0;     % counts number of samples in linear regression piece
numpieces   = 1;     % number of pieces in piecewise linear regression

% go through timeseries, extract voltage peaks & update resistance
nP = 0;              % num voltage/spike peaks we've processed so far
% Create new progress window
waitbar_handles = waitbar(0, [method ' is 0% done'],'Name','rescaleVoltage - close me to stop');
prevProg = 1;
while currt < nT
    % Update Progress window
    try
        waitbar_handles = waitbar(currt/nT, waitbar_handles, [method ' is ', num2str(round( currt/nT*100 )), '% done']);
    catch
        delete(waitbar_handles);
        fprintf(2,['The process has been manually stopped on time: ',num2str(currt*dt),' seconds.\n']);
        break;
    end
    
    % Print current percentage
    if 10*floor(currt/nT*10) ~= prevProg
        fprintf( '\t%s is approx %d%% done\n', method, 10*floor(currt/nT*10) );
        prevProg = 10*floor(currt/nT*10);
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
    % Check three flags: no_glitch, no_noise and not_too_short
    no_glitch = false; 
    no_noise = false;
    not_too_short = false;
    % This has to be checked in a loop cause when one of them happens, the
    % current spike is updated to the next spike, so thresholds and
    % duration have to be checked again.
    while ~(no_glitch && no_noise && not_too_short)
        % if this peak is crazy big assume it's a glitch so move to next spike
        while max( peakfn(vsp) ) > glitchthresh * noisesig
            % move on to next large amplitude event
            startsp = find( peakfn(v( currt:end) ) >= noisethresh, 1, 'first' );
            % no more peaks left in timeseries so exit loop
            if isempty( startsp )
                currt = []; % set to empty since no peaks in remainder of signal
                break;
            end
            prevt   = currt + startsp - 1; % update current time to where next spike is
            [vsp, currt] = getCurrentSpike( peakfn, v, prevt, jumpahead, noisethresh );
            % Make other 2 flags false so it checks for the next chosen
            % spike
            no_noise = false; %#ok<NASGU>
            not_too_short = false; %#ok<NASGU>
        end
        no_glitch = true;

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
            no_glitch = false;
            not_too_short = false; %#ok<NASGU>
        end
        no_noise = true;

        % Check for + and - duration of peak to assess wether it should be
        % considered a spike.
        min_pos_time = 5; % Samples. We can make it a user chosen parameter
        min_neg_time = 5; % Samples. We can make it a user chosen parameter
        sp_too_short = getSpikeDuration(vsp, peakfn, min_pos_time, min_neg_time);
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
            no_noise = false;
            no_glitch = false;
            sp_too_short = getSpikeDuration(vsp, peakfn, min_pos_time, min_neg_time);
        end
        not_too_short = true;
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
    [vmu, vvar] = recursiveUpdate( vmu, vstd^2, lambda, vpeak );
    vstd = sqrt( vvar );
    
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
                % correct( pf, vpeak, Rstd );
                correct( pf, Rtmp, Rstd );
                % Use corrected Rest to update regression coeffs
                Rest  = getStateEstimate(pf);
                
                Rcurr  = Rest(1); % if multiple lagged states get current only
                
                % if not modelling change in resistance as linear
                % regression then force it to change slowly
                % Rest  = Rprev*lambda + Rcurr*(1-lambda);
                
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
                (timeinpiece>=10 && ~iswhite( [Rest_vec(nP-prevpeaks+1:nP-1); Rest_vec(end)] - vpeak_vec(nP-prevpeaks+1:nP), 'bg', 0.05 ) ) ) ...
                || ( Rcurr <= 0 )
            [Rest, Rcoeff, Rmu, Rstd, Rcov] = initRegress(tpeak_vec(nP-prevpeaks+1:nP) - time(currt), ...
                vpeak_vec(nP-prevpeaks+1:nP), lambda, npoles, nzeros);
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
        [Rmu, Rvar] = recursiveUpdate( Rmu, Rvar, lambda, Rcurr );
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
        Rcov_vec(nP,:,:) = Rcov;
        Rmu_vec(nP,1)    = Rmu;
        Rstd_vec(nP,1)   = Rstd;
        vmu_vec(nP,1)    = vmu;
        if ~strcmpi( strip(lower( method )), 'recursive mean' )
            Rcoeff_vec(nP,:) = Rcoeff;
        end
    end
end
% Close progress window
try
    close(waitbar_handles);
catch ME
    if strcmp('Invalid figure handle.',ME.message)
        if exist('waitbar_handles','var')
            clear('waitbar_handles');
        end
    else
        rethrow(ME);
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
        plot(tpeak_vec, noise_vec * glitchthresh, 'r');
        legend('voltage', 'voltage peaks', 'resistance est', 'noise thresh', 'glitch thresh');
        if numel(tpeak_vec) > 1, xlim( [ tpeak_vec(1) tpeak_vec(end) ] );end
    catch ME
        if strcmp('Vectors must be the same length.', ME.message) && isempty(tpeak_vec)
            close;
            error('It looks like there were no spikes found. Check the parameter selection.\nThe value of ''currt'' is: %d', currt);
%             warning('It looks like there were no spikes found. Check the parameter selection.\nThe value of ''currt'' is: %d', currt);
        else
            rethrow(ME);
        end
    end
    
    subplot(222), hold on;
    plot(tpeak_vec, vpeak_vec, 'k.','markersize',marker_size);
    plot(tpeak_vec, Rest_vec,  'm', 'linewidth', 3);
    plot(tpeak_vec, noise_vec * voltoutlier,  'r');
    plot(tpeak_vec, noise_vec * glitchthresh, 'r');
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

% interp doesn't like going outside the input timebounds, so pad Rest
% with end results & pad time to avoid getting NaNs (changes slowly so ok)
if tpeak_vec(1) > time(1)
    tpeak_vec = [time(1); tpeak_vec(:)];
    Rest_vec  = [Rinit; Rest_vec(:)];
end
if tpeak_vec(end) < time(end)
    tpeak_vec = [tpeak_vec(:); time(end)];
    Rest_vec  = [Rest_vec(:); Rest_vec(end)];
end
Rest_vec  = interp1( tpeak_vec, Rest_vec, time );
vrescale  = v(:) ./ Rest_vec(:);

tnow = datetime('now');
str = sprintf("\tFinished rescaling at %s\n", datestr(tnow));
cprintf( 'Keywords', str);

end

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
noisemse    = noisesig;
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
cprintf( 'Keywords', str)


% go through timeseries, extract voltage peaks & update resistance
nP = 0;              % num voltage/spike peaks we've processed so far
while currt < nT
    if mod( nP, 1e3 ) == 0
        fprintf( '\tPS is approx %d%% done\n', round( currt/nT*100 ) );
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
cprintf( 'Keywords', str);


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
jump    = min( remsamp, jumpahead );

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
function too_short = getSpikeDuration(vsp, peakfn, min_pos_time, min_neg_time)
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
    if too_short
        figure; plot(vsp(stidx:endidx));
    end
end

% loop through noise & update mean & mse recursively
function [mu_curr, sig_curr, mse_curr, N_curr] = updateNoiseStats(new_noise, mu_prev, sig_prev, mse_prev, N_prev, lambda)
recursive = false;
N_curr    = N_prev + length( new_noise );

% update noise using either recursive mean or, if assumed stationary,
% using normal mean/sig calculations in real time
if recursive
    % recursive updating of noise, with forgetting
    % init curr vars cuz I seem to have a problem with new noise
    % sometimes being empty?!
    mu_curr = mu_prev; sig_curr = sig_prev; mse_curr = mse_prev;
    for ni=1:length(new_noise)
        % recursive noise estimates
        [mu_curr, var_curr] = recursiveUpdate( mu_prev, sig_prev^2, lambda, new_noise(ni) );
        mse_curr = mse_prev + var_curr;
        mse_curr = mse_prev + (new_noise(ni) - mu_curr) * (new_noise(ni) - mu_prev);
        mu_prev  = mu_curr;
        sig_prev = sqrt(var_curr);
        mse_prev = mse_curr;
    end
    
else
    % normal mean/sig estimates in real time
    mu_curr  = ( mu_prev * N_prev + sum( new_noise ) ) / N_curr;
    mse_curr = mse_prev + sum( (new_noise - mu_curr) .* (new_noise - mu_prev) );
end

sig_curr = sqrt( mse_curr / (N_curr-1) );
end

function [mu_curr, var_curr] = recursiveUpdate( mu_prev, var_prev, lambda, new_obs )
% mu_curr      = mu_prev * (1-dt) + (new_obs - mu_prev) * dt;
mu_curr      = mu_prev  +  (new_obs - mu_prev);
mu_curr      = mu_prev * lambda  +  mu_curr * (1-lambda);
% mu_curr      = mu_prev * lambda  +  (mu_prev + (new_obs - mu_prev)) * (1-lambda);

% S(k) = S(k-1) + (x(k) - M(k-1))*(x(k) - M(k));
% var_curr     = var_prev * dt ...
%             + (new_obs - mu_curr) .* (new_obs - mu_prev) * (1-dt);
var_curr     =  var_prev + ( (new_obs - mu_curr) .* (new_obs - mu_prev) - var_prev );
var_curr     =  var_prev * lambda  +  var_curr * (1-lambda);
% var_curr     = var_prev * lambda ...
%              + ( var_prev + (new_obs - mu_curr) .* (new_obs - mu_prev) ) * (1-lambda);
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
    pf.ResamplingMethod         = 'systematic'; % no idea what this is?!
    %  evolve particles to get new estimates of the resistance
    pf.StateTransitionFcn       = @evolveResistance;
    % obs likelihood - likelihood of each particle state given new obs
    pf.MeasurementLikelihoodFcn = @noiseLikelihood;
end
end

%  initialise with regression on first few samples - else we might miss
%  the first few spikes while the recursive algorithm is converging
function [Rest, Rcoeff, Rmu, Rsig, Rcov] = initRegress(tspike, vspike, lambda, npoles, nzeros, Rcov)
dorecursive = true; % can assume stationary for init section, or do recursively
debug       = true;

% hopefully this will never happen, but just in case!
if length(vspike)==1
    Rcoeff(2) = 0;
    Rcoeff(1) = vspike;
    Rsig = vspike;
    Rmu  = vspike;
    Rcov = eye(2,2) / Rsig;
    Rest = vspike;
    str  = sprintf( 'Init regression for R estimate only has 1 sample!' );
    cprintf( 'Errors', str );
    return;
end

% get init estimate of VAR model to initialise recursive least squares
p      = max(npoles,nzeros);      ndim = npoles + nzeros;
[Rcoeff, ~, vest] = getARcoeff_2vars( vspike, tspike, p );
% NOTE (Artemio): I think this two lines are wrong. I believe a
% is the first 3 elements <p + 1>, and b is the last 2 <p>.
a      =  Rcoeff(2:npoles+1)'; % autoregression coeffs
b      = [Rcoeff(1) Rcoeff(p+2:(p+nzeros))']; % DC + input coeffs
% a = Rcoeff(1:p+1)'; % This is my shot (Artemio)
% b = Rcoeff(p+2:end)'; % This is my shot (Artemio)

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
if nargin<6 || isempty(Rcov)
    Rcov = eye( ndim ) / std(vspike); % inverse cov of predictor matrix
    Rcov = eye( ndim ); % inverse cov of predictor matrix
end
nT     = length(vspike);

Rmu    = mean(vspike); Rmu_prev  = Rmu;
Rvar   = var(vspike);  Rvar_prev = Rvar;
Rest   = vest(end-npoles:end); % init to something just in case insuffic samples
Rsig   = std(vspike);

if debug
    Rmu_vec   = zeros(nT,1);       Rvar_vec  = zeros(nT,1);
    a_vec     = zeros(nT,npoles);  b_vec     = zeros(nT,nzeros);
    Rest_vec  = zeros(nT,1);       Rcov_vec  = zeros(nT,ndim,ndim);
    a_vec(p,:)= a;                 b_vec(p,:)= b;
end
if dorecursive
    for ti=(p+1):nT
        % recursive spike mean/std estimates
        [Rmu, Rvar] = recursiveUpdate( Rmu, Rvar, lambda, vspike(ti) );
        
        % linear regression for spike std
        [a, b, Rcov, Rest] = recursiveLeastSquares( ( tspike(ti-p:ti) - tspike(1) ), ...
            vspike(ti-p:ti), ...
            npoles, nzeros, Rcov, ...
            a, b, 1 );
        
        if debug
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
function [noisemu, noisesig] = initNoiseStdDev(v)
% get initial noise estimate by finding a white time series at the
% beginning of the voltage timeseries
nS   = min(1e4, length(v)); % num init samples to use
vtmp = v(1:nS);
while ~iswhite(vtmp)
    % keep halving number of samples until we get only white samples at
    % the beginning, if we can't find a white subset of samples then just
    % return the first 50
    % Check if 50 is the right number to use. I tried 20 and it didn't work
    nS = ceil( nS/2 );
    if nS < 50
        nS = 50;
        break;
    end
    vtmp = vtmp(1:nS);
end
noisesig =  std( vtmp(1:nS) );
noisemu  = mean( vtmp(1:nS) );
end

function rescaled_voltage = rescaleVoltageVariance(tseries, params)
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