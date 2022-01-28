% [APspikes,APtimes] = extractSpikesUsingTemplates(APtemplates, tseries, method, params)
% 
% Extract spikes that match the input AP templates. For each template we
% match on shape, but then within a template we allow different families
% based on the amplitude of the spikes. That way multiple axons can have
% the same shape, but hopefully can be separated out using amplitude. 
% Inputs:
%  APtemplates - array of spike shapes for the templates, all of same length
%                format: template_length  x  num_templates
%  tseries     - tseries struct from SET gui of type voltage, to find spikes in 
%  method      - method used to extract potential spikes
%                'matched filter', or 'wavelets'
%  params      - parameters for extracting spikes, set using SET gui
% Outputs:
%  APspikes    - spike for each AP template, organised by family (scale)
%                within each template, so that there's a matrix of 
%                time x num_families within each template cell
%  APtimes     - the time of each peak for the spikes within a family,
%                within a template, so that APtimes has a cell for each AP
%                template, and within the cell another cell for each
%                family, and within that cell an array of peak times for
%                each spike within that family

% TO DO: if spike peaks aren't white then it suggests an unmodelled change,
% such as what I'm seeing which is sudden jumps in amplitude that are then
% constant for a while, then another jump down in amplitude
function [APspikes, APtimes] = extractSpikesUsingTemplates( APtemplates, APnumsamples, tseries, method, params, normAPs ,opts)
   
   if exist('opts', 'var') && isfield(opts,'loadingWindowOn'), progress_window = opts.loadingWindowOn; 
   else, progress_window = false; end
   if opts.auto_params
      % Get auto params
      pp = getAutomaticParams(tseries, [], 'extractspikes', method, params);
      % If returned value is empty, leave the default parameters untouched
      if ~isempty(pp)
         params = pp;
      end
   end  
   % AP templates can have different lengths if in a cell, else all forced
   % to same length if in a matrix
   if iscell( APtemplates ) 
      uniqueAPLength = true; 
      nAP       = length(APtemplates); % num samples in each template, & num templates
      nT        = cellfun( @(t) size( t, 1 ), APtemplates );
      peakind   = cellfun( @(t) getMaxInd( t, 1 ), APtemplates );
      peakN     = round( median( peakind ) ); % mean index of template peaks
      % peak function calculates the peak of the spike - e.g. using total diff
      % btwn min & max, or just max
      peakfn    = @(sp) max(sp) - min(sp); 
      peakfn    = @(sp) max(sp); 
      peakfn    = @(sp) [ min(sp) max(sp) ]; 
   else
      uniqueAPLength = false; 
      [nT, nAP] = size(APtemplates); % num samples in each template, & num templates
      peakind   = getMaxInd( APtemplates, 1 );
      peakN     = round( median( peakind ) ); % mean index of template peaks
      % peak function calculates the peak of the spike - e.g. using total diff
      % btwn min & max, or just max
      peakfn    = @(sp) sp(peakN) - min(sp); 
      peakfn    = @(sp) sp(peakN); 
      peakfn    = @(sp) [ min(sp) sp(peakN) ]; 
   end

   APspikes  = [];  APtimes = [];
   orig_nAP  = nAP; % if allowing new templates, record how many we started with
   
   % peakdiff function calculates diff btwn peaks, e.g. btwn total range or
   % sums diff btwn maxes & mins, & perhaps accounts for time btwn spikes
   % and expects inputs of the form [peak_neg  peak_pos]
   % - weight positive peak more heavily than minimum peak
	% peakdiff  = @(curr_peak, last_peaks) sum( abs((curr_peak - last_peaks) .* [1/2 2] ./ last_peaks), 2 );
	% peakdiff  = @(x, y) abs(((x(:,2)-x(:,1))-(y(:,2)-y(:,1)))./(y(:,2)-y(:,1))) ;

   % divide by 2 to avg diff across pos & neg peaks
   peakdiff  = @(curr_peak, last_peaks) sum( abs((curr_peak - last_peaks) .* [1 1] ./ last_peaks), 2 ) / 2;
   timefn    = @(sp, st) st( getMaxInd(sp) );
   
   switch lower(method)
      case 'matched filter'
         if uniqueAPLength
            [spikes, stimes, spikesFull, stimesFull] = getSpikesByThresholding( tseries, params );
            tmpsp  = spikes;
            tmpt   = stimes;
            spikes = spikesFull;
            stimes = stimesFull;
         else
            [spikes, stimes, spikesFull] = getSpikesByThresholding( tseries, params, nT, peakN );
         end

      % AP templates as mother wavelets gave terrible results, & since we
      % need spikes per template, ditch the use of wavelets for extracting
      % spikes (you can still use them to generate the templates though)
%       case 'wavelets'
%          str = sprintf('Not yet, will come in time...');
%          displayErrorMsg(str);
%          return;
%          
      otherwise
         str = sprintf('What''s this now? Don''t know what you''re talking about...');
         displayErrorMsg(str);
         return;
   end
%%
   % now that we have a list of possible spikes matching user's input
   % constraints, we need to sort the spikes into families. Use the AP
   % templates to match shape, but sort into separate families according to
   % scale, so that each template can have multiple families
   matchtype   = params.match_type.value;
   matchthresh = params.match_similarity.value;
   APfamilies  = cell(nAP,1); % cell for each template, with a cell for each family inside
   newTemplates= params.allow_new_aps.value;      % allow new AP templates or not
   lambda      = params.forgetting_factor.value;
   kappa_pos   = params.kappa_pos.value; % Times variance
   kappa_neg   = params.kappa_neg.value; % Times variance
      
   % Print progress
   str = sprintf( '\tSorting spikes...\n' );
   cprintf( 'Keywords', str );   
   nS  = size(spikes, 2);
   prevProg = 1;
   % Create new progress window
   if progress_window
      waitbar_handles = waitbar(0,[method ' is 0% done'],...
         'Name','Extracting spikes - close me to stop');
   end
   
   switch matchtype
         case 'corr'
            matchfn = @( sp, temp ) corr( temp, sp );
            matchfn = @( sp, temp ) max( abs( xcorr( temp, sp ) ) );

         case 'cov'
            % to compare using covariance means that the APs have to have
            % the same amplitude to be considered part of the same family.
            % In order to have a match threshold independent of
            % amplitude, scale by the size of the template's signal. 
            % Since cov is voltage squared, scale by variance rather than 
            % std dev. Scale by var of templates, and of current spike,
            % and then take the smallest result. This means that size is
            % taken care of, and avoids the situation where cov is
            % approx 1 but the fit is bad, but just happened to scale to
            % 1 because of the sum in the numerator
            matchfn = @( sp, temp ) min( cov_x_Yvec( sp, temp ) ./ var( temp ), cov_x_Yvec( sp, temp ) ./ var( sp ) ); 

      otherwise
            str = sprintf('%s is an unknown match type', matchtype);
            displayErrorMsg(str);
            return;
   end
   
   try
      for si=1:nS % for each possible spike (format: time x spike)
         % Print current percentage
         if 10*floor(si/nS*10) ~= prevProg
             fprintf( '\t%s is approx %d%% done\n', method, 10*floor(si/nS*10) );
             prevProg = 10*floor(si/nS*10);
         end

         % Update Progress window
         if progress_window
            try
               waitbar_handles = waitbar(si/nS, waitbar_handles,...
                  [method ' is ', num2str(round( si/nS*100 )), '% done']);
            catch
               delete(waitbar_handles);
               str = ['The process has been manually stopped on time: ',num2str(stimes{si}),' seconds.\n'];
               printMessage('off','Error',str);
               break;
            end
         end

         if uniqueAPLength 
            curr_spike = spikes{si};
         else
            curr_spike = spikes(:,si);
         end
         curr_stime = stimes{si};
         rho = spikeCorrWithTemplates( curr_spike, APtemplates, matchfn );

         % If spike is sufficiently similar in shape to a template, 
         % add it to one of the template families, or start a new family 
         % if the spike amplitude is not sufficiently similar.
         % - gotta account for cov so rho has to be within band around 1
         if any( rho > matchthresh )
            k = compare_templates(rho, curr_spike, APtemplates, APfamilies); % template with closest match

            % if no spike families yet, create one
            if isempty(APfamilies{k})
               % Distribution of peaks (mean and var)
               [fam_mu, fam_var] = initialize_mu( spikes, si, matchthresh, peakfn, peakdiff ); 
               APfamilies{k}{1}  = newAxonFamily( curr_spike, peakfn, curr_stime, timefn, fam_mu, fam_var, opts );

            % if we have spike families see if we're at the right
            % scale, so get size of current peak, & time diff btwn
            % last spike and this spike, & check if the peak is
            % sufficiently the same to be considered part of the fam
            else
               curr_peak = peakfn(curr_spike);           
               varn      = getStructFieldFromCell( APfamilies{k}, 'varn', 'array' );
               varp      = getStructFieldFromCell( APfamilies{k}, 'varp', 'array' );
               mun       = getStructFieldFromCell( APfamilies{k}, 'mun', 'array' );
               mup       = getStructFieldFromCell( APfamilies{k}, 'mup', 'array' );

               % get rate of change from last peak (this is populated with the
               % mean spike now) to current peak, but then switch around
               % because I was finding that spikes were always shrinking 
               % over time, because if last_peaks is smaller (as in we're 
               % growing) when we divide by last peaks you're dividing by a 
               % small number, which makes the change rate larger, & less
               % likely to be within the allowed limit
               change_rates1 = peakdiff( curr_peak, [mun mup] );
               change_rates2 = peakdiff( [mun mup], curr_peak );

               % use the change_rates with the smallest value since it'll be
               % the min change rate that's selected below
               change_rates  = ternaryOp( min( change_rates1(:) ) < min( change_rates2(:) ), ...
                                          change_rates1, change_rates2 );

               % spike peak (pos & neg) must be within range of a gaussian
               % distribution's mu and variance
               stdn_  = sqrt(varn);
               stdp_  = sqrt(varp);
               valid1 = zeros(size(mup));
               valid2 = zeros(size(mun));
               for ti = 1:numel(mup)
                  valid1(ti) = (curr_peak(1) >= mun(ti) - kappa_neg * stdn_(ti)) && (curr_peak(1) <= mun(ti) + kappa_neg * stdn_(ti));
                  valid2(ti) = (curr_peak(2) >= mup(ti) - kappa_pos * stdp_(ti)) && (curr_peak(2) <= mup(ti) + kappa_pos * stdp_(ti));
               end
               valid = valid1 & valid2;

               if any(valid)
                  % Adding one spike to axon family of best fit
                  [minrate, mini] = min( change_rates(valid) );
                  vind = find( valid );
                  fi   = vind( mini );
                  fam  = APfamilies{k}{fi};

                  % if current spike peak is within required rate of
                  % change in peak amplitude for a spike family, add it
                  % to the fam. Else start a new family.              
                  if uniqueAPLength
                     [ fam.meanspike, fam.spikes, curr_stime, fam.stimes ] = addSpikeToTemplate( fam.meanspike, curr_spike, fam.spikes, curr_stime, fam.stimes );
                  else
                     fam.spikes(:,end+1) = curr_spike;
                     fam.stimes(:,end+1) = curr_stime;
                     N = size( fam.spikes, 2 ); % update number of spikes
                     % efficient calculation of meanspike for large families
                     fam.meanspike = ( fam.meanspike * (N-1) + curr_spike ) / N;
                     % APfamilies{k}{i}.meanspike = APfamilies{k}{end}.meanspike + curr_spike * (1-lambda);
                  end

                  % Distribution of peaks (mean and var)
                  [mup_curr, varp_curr] = recursiveUpdate( fam.mup, fam.varp, curr_peak(2), lambda );
                  [mun_curr, varn_curr] = recursiveUpdate( fam.mun, fam.varn, curr_peak(1), lambda );
                  fam.mun       = mun_curr;
                  fam.varp      = varp_curr;
                  fam.varn      = varn_curr;             
                  fam.last_peak = peakfn( curr_spike ); % size of last spike in family
                  fam.last_time = timefn( curr_spike, curr_stime ); % time of last spike             
                  if ~strcmp( 'none', opts.debugOption )
                     fam   = updateDebugInfo( fam, mup_curr, mun_curr, varp_curr, varn_curr );
                  end
                  APfamilies{k}{fi} = fam;
               else
                  % Creating new axon family
                  [fam_mu, fam_var]     = initialize_mu( spikes, si, matchthresh, peakfn, peakdiff ); 
                  APfamilies{k}{end+1}  = newAxonFamily( curr_spike, peakfn, curr_stime, timefn, fam_mu, fam_var, opts );
               end
            end

            % don't update template because we've probably already included
            % these spikes, unless the template is a new one
            if k > orig_nAP
               % update template by incorporating current spike into
               % mean - we don't care about amplitue
               if normAPs
                  curr_spike  = curr_spike / std(curr_spike); 
               end
               % Equal importance to all past spikes
               % APtemplates(:,k) = (APtemplates(:,k)*APnumsamples(k) + ...
               % curr_spike) / (APnumsamples(k)+1); % uncomment this line and
               % comment the next one to go back
               % Forgetting factor - lambda defines how much weight do older spikes
               % have
               APtemplates(:,k) = APtemplates(:,k) * lambda  +  curr_spike * (1 - lambda);
               APnumsamples(k)  = APnumsamples(k) + 1;
            end

         % Spike is not sufficiently similar to one of the families
         % so start a new AP template if user has allowed new ones
         else
            if newTemplates
               % make a new AP family within the templates
               if uniqueAPLength 
                  APtemplates{end+1} = curr_spike;
               else
                  APtemplates(:, end+1) = curr_spike;
               end
               APnumsamples(end+1)   = 1;
               % Creating new axon family
               [fam_mu, fam_var]     = initialize_mu( spikes, si, matchthresh, peakfn, peakdiff ); 
               APfamilies{end+1}{1}  = newAxonFamily( curr_spike, peakfn, curr_stime, timefn, fam_mu, fam_var, opts );
               nAP = nAP + 1;
            end
         end
      end
   catch ME
      str = getCatchMEstring( ME, 'main: ' );
      cprintf( 'Keywords*', str );
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

   try
      % get rid of families with very few spikes
      for fi=1:length( APfamilies )
         if ~isempty(APfamilies{fi})
             nf = cellfun( @(fam) size( fam.spikes,2 ), APfamilies{fi} );
             toosmall = nf <= params.min_spiking_threshold.value; 
             APfamilies{fi}( toosmall ) = [];
         end
      end
      % get rid of empty families, just in case they snuck in there!
      nf = cellfun( @length, APfamilies );
      loner = nf <= 0; % get empty families
      APtmp = APfamilies; % make a copy before deleting small families
      APfamilies( loner ) = [];
      if isempty( APfamilies )
         displayErrorMsg( 'No axon families found with this parameter configuration' );
         return; 
      end
      nAP = length( APfamilies ); % update num families after removing loners

      if ~strcmp( 'none', opts.debugOption )
         figure; 
         % if not too many templates then plot each templates mean & var in
         % a different column
         if nAP<4
            nc = nAP; nr = 2;
         else
            % too many templates --> just plot in pairs
            nc = ceil( sqrt( nAP*2 ) ); nc = ternaryOp( round(nc/2)~=nc/2, nc+1, nc );
            nr = ceil( nAP*2 / nc );
         end
         for ti = 1:size(APfamilies) % for each template
            mu_pos = -Inf; var_pos = -Inf;
            mu_neg =  Inf; var_neg = -Inf;
            nfam   = size( APfamilies{ti}, 2 );
            cols   = getColourMatrix( nfam );
            for fi = 1:nfam % for each family within the template
               subplot(nr,nc,ti); hold on;
                  plot( APfamilies{ti}{fi}.stimes(1,:), APfamilies{ti}{fi}.mup_vector, 'color', cols(fi,:) );
                  plot( APfamilies{ti}{fi}.stimes(1,:), APfamilies{ti}{fi}.mun_vector, 'color', cols(fi,:) );
                  mu_pos = ternaryOp( mu_pos < max( APfamilies{ti}{fi}.mup_vector ), max( APfamilies{ti}{fi}.mup_vector ), mu_pos );
                  mu_neg = ternaryOp( mu_neg > min( APfamilies{ti}{fi}.mun_vector ), min( APfamilies{ti}{fi}.mun_vector ), mu_neg );
               subplot(nr,nc,nAP+ti); hold on;
                  plot( APfamilies{ti}{fi}.stimes(1,:), APfamilies{ti}{fi}.varp_vector, 'color', cols(fi,:), 'linewidth', 1 );
                  plot( APfamilies{ti}{fi}.stimes(1,:), APfamilies{ti}{fi}.varn_vector, '--', 'color', cols(fi,:), 'linewidth', 1 );
                  var_pos = ternaryOp( var_pos < max( APfamilies{ti}{fi}.varp_vector ), max( APfamilies{ti}{fi}.varp_vector ), var_pos );
                  var_neg = ternaryOp( var_neg < max( APfamilies{ti}{fi}.varn_vector ), max( APfamilies{ti}{fi}.varn_vector ), var_neg );
            end
            subplot(nr,nc,ti),     ylim( [ mu_neg,  mu_pos ] );
            subplot(nr,nc,nAP+ti), ylim( [      0, var_pos ] );
         end
      end


      % Each valid spike has been allocated to an AP template, & to an AP
      % familiy within that. Now for each AP family we need to construct a
      % timeseries
      if ~params.plot_spikes_consecutively.value
         for ti=1:nAP % for each AP template
            fam = APfamilies{ti}; 
            nF  = length(fam);

            try
               % try allocating a full tseries for each family in template
               fam_tseries = zeros(length(tseries.data),nF);
               fam_stimes  = cell(nF, 1);

               for ff=1:nF % for each family within the AP templates
                  % if we're creating the timeseries extract spikes at
                  % appropriate times & copy into voltage timeseries
                  sind      = round(fam{ff}.stimes(:) / tseries.dt);
                  valid     = ~isnan( sind );
                  invalid   = isnan( sind );
                  fam_tseries(valid,ff) = fam{ff}.spikes(valid);
                  peakind   = round( median( getMaxInd( fam{ff}.spikes ) ) );
                  fam_stimes{ff}  = fam{ff}.stimes(peakind,:);
               end
               APspikes{ti} = fam_tseries;
               APtimes{ti}  = fam_stimes;

            catch ME
               % running outta memory so treat each family separately so it's
               % event based rather than disrete time
                 % contains( ME.identifier, 'sizelimit', 'IgnoreCase', true)
               if checkThrownError( ME, 'SizeLimitExceeded' )
                  for ff=1:nF 
                     % for each family within the AP templates, extract the times of
                     % spike peaks
                     peakind = round( median( getMaxInd( fam{ff}.spikes ) ) );
                     fam_stimes{ff}  = fam{ff}.stimes(peakind,:);
                  end
                  fam{1}.time     = tseries.time;
                  APspikes{ti} = fam;
                  APtimes{ti}  = fam_stimes;
                  
               else
                  str = getCatchMEstring( ME, 'extractSpikesUsingTemplates: error generating spike timeseries\n\t', false );
                  runtimeErrorHandler(ME);
               end
            end
         end

      % squish spikes together so they can be seen as one long series of spikes
      else
         % determine which family's spikes take longest to run 
         % consecutively because we'll make all timeseries accommodate 
         % that length (currently time limits are common to whole dataset)
         % max duration of consecutive spikes, i.e. largest time x num_spikes
         maxtime   = max( cellfun(@(template) max( cellfun(@(family) numel( family.stimes ), template ) ), APfamilies ) );
         for ti=1:nAP % for each AP template
            fam = APfamilies{ti}; 
            nF  = length(fam);
            fam_tseries = zeros( maxtime, nF ); % time_x_spikes x num families
            fam_stimes  = cell(nF, 1);
            for ff=1:nF % for each family within the AP templates
               [ nT, nSp ] = size( fam{ff}.stimes ); % number of spikes for this family
               % populate family timeseries by running spikes one after the other
               fam_tseries(1:(nSp*nT), ff) = toVec(fam{ff}.spikes(:));
               peakind = floor( median( getMaxInd( fam{ff}.spikes ) ) );
               fam_stimes{ff} = fam{ff}.stimes(peakind,:);
            end
            APspikes{ti} = fam_tseries;
            APtimes{ti}  = fam_stimes;
         end

      end

   catch ME
      str = getCatchMEstring( ME, 'end: ' );
   end
end

function [mu, sig] = initialize_mu( spikes, si, matchthresh, peakfn, diffpeakfn )
   maxlag = 2; % number of lags in cross-corr
   % Initialize the mean peak for a new family, triggered by spike at index 
   % si not having a family of similar size. Goes through the first N
   % spikes that match the template and chooses that value as the mean.
   if iscell( spikes ), uniqueAPLength = true; else, uniqueAPLength = false; end
   if uniqueAPLength
      nSp = length( spikes );
      fsp = spikes{si}; % spike we're making the family from
   else
      nSp = size( spikes, 2 );
      fsp = spikes(:,si);
   end
   fam_peaks = [ min(fsp) max(fsp) ];
   si  = si + 1; % start matching with spikes after current spike
   allowed_diff = 0.5; % Allowed difference between spikes to initialize mu
   % if we've run outta spikes in the spike train, return
   if si >= nSp
      mu = peakfn( fsp ); sig = 0; 
      return; 
   end
   
   cnt = 0;
   pk = [];
   while cnt < 10 && si < nSp
      if uniqueAPLength
         sp = spikes{si}; % spike we're making the family from
      else
         sp = spikes(:,si);
      end
      sp_peaks = [ min(sp) max(sp) ];
      % align the spikes at their peaks 
      [ sp, fsp ] = alignDiffLengthSpikes( sp, fsp );
      c   = xcorr( sp, fsp, maxlag, 'coeff' );
      rho = max( abs(c) );
      d1  = diffpeakfn( sp_peaks, fam_peaks ); % diff in peaks btwn family spike & sp
      d2  = diffpeakfn( fam_peaks, sp_peaks ); % diff in peaks in reverse order
      d   = min( d1, d2 ); % get largest diff btwn pos & neg peak changes
%       d   = max( d1, d2 ); % get largest diff btwn pos & neg peak changes

      % If spike is closeley correlated to the template
      if rho > matchthresh && d < max(fsp)*0.3 % TO DO: DON'T HARD-CODE ALLOWABLE DIFFERENCE !! 
         cnt     = cnt + 1;
         if uniqueAPLength
            pk_  = peakfn( spikes{si} );
         else
            pk_  = peakfn(spikes(:,si));
         end
         pk(cnt,:) = pk_;
      end
      si = si + 1;
   end
   
   % If at least 2 spikes matched the template (hopefully 5+ did!)
   if cnt > 2
      mu  = mean(pk, 1);
      sig = var(pk, [], 1 );
   else
      fsp = spikes(:,si);
      mu  = peakfn( fsp );
      sig = sqrt( max( fsp ) - min( fsp ) ) * ones(1,2);
   end
end

% Spawn a new axon family with all the necessary struct fields
function family = newAxonFamily( curr_spike, peakfn, curr_stime, timefn, fam_mu, fam_var, opts )
   family.spikes    = curr_spike;
   family.stimes    = curr_stime;
   family.meanspike = curr_spike;
   family.last_peak = peakfn( curr_spike ); % size of last spike in family
   family.last_time = timefn( curr_spike, curr_stime ); % time of last spike
   family.mun       = fam_mu(1); 
   family.mup       = fam_mu(2); 
   family.varn      = fam_var(1);
   family.varp      = fam_var(2);
  
   if ~strcmp( 'none', opts.debugOption )
      family.mup_vector  = fam_mu(2);
      family.mun_vector  = fam_mu(1);
      family.varp_vector = fam_var(2);
      family.varn_vector = fam_var(1);
   end
end

function family = updateDebugInfo( family, mup_curr, mun_curr, varp_curr, varn_curr )
   family.mup_vector  = [ family.mup_vector mup_curr ];
   family.mun_vector  = [ family.mun_vector mun_curr ];
   family.varp_vector = [ family.varp_vector sqrt( varp_curr ) ];
   family.varn_vector = [ family.varn_vector sqrt( varn_curr ) ];
end



