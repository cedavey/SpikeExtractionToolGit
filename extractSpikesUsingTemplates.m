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
function [APspikes, APtimes] = extractSpikesUsingTemplates( APtemplates, APnumsamples, tseries, method, params, normAPs ,opts)
   
   if isfield(opts,'loadingWindowOn'), progress_window = opts.loadingWindowOn; else, progress_window = false; end
   if opts.auto_params
      % Get auto params
      pp = getAutomaticParams(tseries, [], 'extractspikes', method, params);
      % If returned value is empty, leave the default parameters untouched
      if ~isempty(pp)
         params = pp;
      end
   end   

   APspikes  = [];  APtimes = [];
   [nT, nAP] = size(APtemplates); % num samples in each template, & num templates
   orig_nAP  = nAP; % if allowing new templates, record how many we started with
   peakind   = getMaxInd( APtemplates, 1 );
   peakN     = round( mean( peakind ) ); % mean index of template peaks
   % peak function calculates the peak of the spike - e.g. using total diff
   % btwn min & max, or just max
   peakfn    = @(sp) sp(peakN) - min(sp); 
   peakfn    = @(sp) sp(peakN); 
   peakfn    = @(sp) [ min(sp) sp(peakN) ]; 
   
   % peakdiff function calculates diff btwn peaks, e.g. btwn total range or
   % sums diff btwn maxes & mins, & perhaps accounts for time btwn spikes
   % and expects inputs of the form [peak_neg  peak_pos]
   % - weight positive peak more heavily than minimum peak
	% peakdiff  = @(curr_peak, last_peaks) sum( abs((curr_peak - last_peaks) .* [1/2 2] ./ last_peaks), 2 );
	% peakdiff  = @(x, y) abs(((x(:,2)-x(:,1))-(y(:,2)-y(:,1)))./(y(:,2)-y(:,1))) ;

   % divide by 2 to avg diff across pos & neg peaks
   peakdiff  = @(curr_peak, last_peaks) sum( abs((curr_peak - last_peaks) .* [1 1] ./ last_peaks), 2 ) / 2;
   timefn    = @(st) st(peakN);
   
   switch lower(method)
      case 'matched filter'
         [spikes, stimes] = getSpikesByThresholding(tseries, params, nT, peakN);

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
      
      curr_spike = spikes(:,si);
      curr_stime = stimes{si};
      switch matchtype
            case 'corr'
               % rho = corr(APtemplates, curr_spike); % This line was not
               % being used. I commented it to save time. It was taking
               % several minutes to run it in every iteration.
               
               % use autocorrelation to allow spikes to not be perfectly
               % aligned
               rho = zeros( nAP, 1 );
               for ap=1:nAP
                  % Replaced the option 'normalized' with 'coeff' to ensure
                  % compatibility with older versions. New versions accept
                  % both.
                  [c, lags] = xcorr( curr_spike, APtemplates(:,ap), 2, 'coeff');
                  % use max autocorr
                  rho(ap) = max( abs(c) );
               end
               
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
               rho = min( [ cov_x_Yvec(curr_spike, APtemplates) ./ var(APtemplates); cov_x_Yvec(curr_spike, APtemplates) ./ var(curr_spike) ] );
               
         otherwise
               str = sprintf('%s is an unknown match type', matchtype);
               displayErrorMsg(str);
               return;
      end
      
      % If spike is sufficiently similar in shape to a template, 
      % add it to one of the template families, or start a new family 
      % if the spike amplitude is not sufficiently similar.
      % - gotta account for cov so rho has to be within band around 1
      if any( rho > matchthresh )
         k = compare_templates(rho, curr_spike, APtemplates, APfamilies); % template with closest match
         
         % if no spike families yet, create one
         if isempty(APfamilies{k})
            APfamilies{k}{1}.spikes    = curr_spike;
            APfamilies{k}{1}.stimes    = curr_stime;
            APfamilies{k}{1}.meanspike = curr_spike;
            APfamilies{k}{1}.last_peak = peakfn(curr_spike); % size of last spike in family
            APfamilies{k}{1}.last_time = timefn(curr_stime); % time of last spike
            % Distribution of peaks (mean and var)
            [fam_mu, fam_var] = initialize_mu( spikes, si, matchthresh, peakfn, peakdiff ); 
             APfamilies{k}{1}.mun  = fam_mu(1); 
             APfamilies{k}{1}.mup  = fam_mu(2); 
             APfamilies{k}{1}.varn = fam_var(1);
             APfamilies{k}{1}.varp = fam_var(2);

            if ~strcmp('none',opts.debugOption)
               APfamilies{k}{1}.mup_vector  = [APfamilies{k}{1}.mup];
               APfamilies{k}{1}.mun_vector  = [APfamilies{k}{1}.mun];
               APfamilies{k}{1}.varp_vector = [APfamilies{k}{1}.varp];
               APfamilies{k}{1}.varn_vector = [APfamilies{k}{1}.varn];
            end
            
         % if we have spike families see if we're at the right
         % scale, so get size of current peak, & time diff btwn
         % last spike and this spike, & check if the peak is
         % sufficiently the same to be considered part of the fam
         else
            curr_peak     = peakfn(curr_spike);
            
            varn          = getStructFieldFromCell(APfamilies{k}, 'varn', 'array');
            varp          = getStructFieldFromCell(APfamilies{k}, 'varp', 'array');
            mun           = getStructFieldFromCell(APfamilies{k}, 'mun', 'array');
            mup           = getStructFieldFromCell(APfamilies{k}, 'mup', 'array');

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
            for j = 1:numel(mup)
               valid1(j) = (curr_peak(1) >= mun(j) - kappa_neg * stdn_(j)) && (curr_peak(1) <= mun(j) + kappa_neg * stdn_(j));
               valid2(j) = (curr_peak(2) >= mup(j) - kappa_pos * stdp_(j)) && (curr_peak(2) <= mup(j) + kappa_pos * stdp_(j));
            end
            valid = valid1 & valid2;
                        
            if any(valid)
               % Adding one spike to current axon family              
               [minrate, mini] = min( change_rates(valid) );
               % if current spike peak is within required rate of
               % change in peak amplitude for a spike family, add it
               % to the fam. Else start a new family.
               valid_ind = find( valid );
               i = valid_ind( mini );
               APfamilies{k}{i}.spikes(:,end+1) = curr_spike;
               APfamilies{k}{i}.stimes(:,end+1) = curr_stime;
               N = size( APfamilies{k}{i}.spikes, 2 ); % update number of spikes
               % efficient calculation of meanspike for large families
               APfamilies{k}{i}.meanspike = ( APfamilies{k}{i}.meanspike * (N-1) + curr_spike ) / N;
               % APfamilies{k}{i}.meanspike = APfamilies{k}{end}.meanspike + curr_spike * (1-lambda);

               % Distribution of peaks (mean and var)
               mup_prev  = APfamilies{k}{i}.mup;
               mun_prev  = APfamilies{k}{i}.mun;
               varp_prev = APfamilies{k}{i}.varp;
               varn_prev = APfamilies{k}{i}.varn;
               [mup_curr, varp_curr] = recursiveUpdate( mup_prev, varp_prev, curr_peak(2), lambda );
               [mun_curr, varn_curr] = recursiveUpdate( mun_prev, varn_prev, curr_peak(1), lambda );
               APfamilies{k}{i}.mun  = mun_curr;
               APfamilies{k}{i}.varp = varp_curr;
               APfamilies{k}{i}.varn = varn_curr;
               
               if ~strcmp( 'none', opts.debugOption )
                  APfamilies{k}{i}.mup_vector  = [ APfamilies{k}{i}.mup_vector mup_curr];
                  APfamilies{k}{i}.mun_vector  = [ APfamilies{k}{i}.mun_vector mun_curr];
                  APfamilies{k}{i}.varp_vector = [ APfamilies{k}{i}.varp_vector sqrt(APfamilies{k}{i}.varp)];
                  APfamilies{k}{i}.varn_vector = [ APfamilies{k}{i}.varn_vector sqrt(APfamilies{k}{i}.varn)];
               end
               
               % getting peak to peak from family mean now, because it
               % wasn't working at all well using the last spike
               APfamilies{k}{i}.last_peak   = peakfn( curr_spike ); % size of last spike in family
               APfamilies{k}{i}.last_time   = timefn( curr_stime ); % time of last spike
               
            else
               % Creating new axon family
               APfamilies{k}(end+1)         = cell(1);
               APfamilies{k}{end}.spikes    = curr_spike;
               APfamilies{k}{end}.stimes    = curr_stime;               
               APfamilies{k}{end}.meanspike = curr_spike; 
               APfamilies{k}{end}.last_peak = peakfn( curr_spike ); % size of last spike in family
               APfamilies{k}{end}.last_time = timefn( curr_stime ); % time of last spike
                
               [fam_mu, fam_var] = initialize_mu( spikes, si, matchthresh, peakfn, peakdiff ); 
                APfamilies{k}{end}.mun  = fam_mu(1); 
                APfamilies{k}{end}.mup  = fam_mu(2); 
                APfamilies{k}{end}.varn = fam_var(1);
                APfamilies{k}{end}.varp = fam_var(2);
                              
               if ~strcmp( 'none', opts.debugOption )
                  APfamilies{k}{end}.mup_vector  = fam_mu(2);
                  APfamilies{k}{end}.mun_vector  = fam_mu(1);
                  APfamilies{k}{end}.varp_vector = sqrt( fam_var(2) );
                  APfamilies{k}{end}.varn_vector = sqrt( fam_var(1) );
               end
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
            APtemplates(:, end+1) = curr_spike;
            APnumsamples(end+1)   = 1;
            APfamilies(end+1,1)   = cell(1);
            APfamilies{end}{1}.spikes    = curr_spike;
            APfamilies{end}{1}.stimes    = curr_stime;
            APfamilies{end}{1}.last_peak = peakfn(curr_spike); % size of last spike in family
            APfamilies{end}{1}.last_time = timefn(curr_stime); % time of last spike
            nAP = nAP + 1;
         end
      end
   end
catch ME
   str = getCatchMEstring( ME, 'main: ' );
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
   if ~strcmp('none',opts.debugOption)
      figure;subplot(1,2,1); hold;
      for j = 1:size(APfamilies)
         for jj = 1:size(APfamilies{j},2)
            plot(APfamilies{j}{jj}.stimes(1,:),APfamilies{j}{jj}.mup_vector);
         end
      end
      subplot(1,2,2); hold
      for j = 1:size(APfamilies)
         for jj = 1:size(APfamilies{j},2)
            plot(APfamilies{j}{jj}.stimes(1,:),APfamilies{j}{jj}.varp_vector,'--');
         end
      end
   end
   
   % get rid of families with very few spikes
   for fi=1:length( APfamilies )
      if ~isempty(APfamilies{fi})
          nf = cellfun( @(fam) size( fam.spikes,2 ), APfamilies{fi} );
          toosmall = nf <= params.min_spiking_threshold.value; 
          APfamilies{fi}(toosmall) = [];
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
               sind = round(fam{ff}.stimes(:) / tseries.dt);
               fam_tseries(sind,ff) = fam{ff}.spikes(:);
               fam_stimes{ff}  = fam{ff}.stimes(peakN,:);
            end
            APspikes{ti} = fam_tseries;
            APtimes{ti}  = fam_stimes;
            
         catch ME
            % running outta memory so treat each family separately so it's
            % event based rather than disrete time
            for ff=1:nF 
               % for each family within the AP templates, extract the times of
               % spike peaks
               fam_stimes{ff}  = fam{ff}.stimes(peakN,:);
            end
            fam{1}.time     = tseries.time;
            APspikes{ti} = fam;
            APtimes{ti}  = fam_stimes;
         end
      end
      
   % squish spikes together so they can be seen as one long series of spikes
   else
      % determine which family has the most spikes, because we'll make all
      % timeseries accommodate that length (currently time limits are
      % common to whole dataset)
      maxspikes = max( cellfun(@(template) max( cellfun(@(family) length(family.stimes), template ) ), APfamilies ) );
      for ti=1:nAP % for each AP template
         fam = APfamilies{ti}; 
         nF  = length(fam);
         fam_tseries = zeros(maxspikes * nT, nF); % max spikes*spike length x num families
         fam_stimes  = cell(nF, 1);
         for ff=1:nF % for each family within the AP templates
            nspikes = size(fam{ff}.stimes, 2); % number of spikes for this family
            % populate family timeseries by running spikes one after the other
            fam_tseries(1:(nspikes*nT), ff) = toVec(fam{ff}.spikes(:));
            fam_stimes{ff} = fam{ff}.stimes(peakN,:);
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
   % spikes that match the template and choses that value as the mean.
   fsp = spikes(:,si); % spike we're making the family from
   fam_peaks = [ min(fsp) max(fsp) ];
   si  = si + 1; % start matching with spikes after current spike
   % if we've run outta spikes in the spike train, return
   if si >= size( spikes, 2 )
      mu = peakfn( fsp ); sig = 0; 
      return; 
   end
   
   cnt = 0;
   pk = [];
   while cnt < 10 && si < size(spikes,2)
      sp  = spikes(:,si);
      sp_peaks = [ min(sp) max(sp) ];
      c   = xcorr( sp, fsp, maxlag, 'coeff' );
      rho = max( abs(c) );
      d1  = diffpeakfn( sp_peaks, fam_peaks ); % diff in peaks btwn family spike & sp
      d2  = diffpeakfn( fam_peaks, sp_peaks ); % diff in peaks in reverse order
      d   = min( d1, d2 ); % get largest diff btwn pos & neg peak changes

      % If spike is closeley correlated to the template
      if rho > matchthresh && d < 0.3 % TO DO: DON'T HARD-CODE ALLOWABLE DIFFERENCE !! 
         cnt     = cnt + 1;
         pk_     = peakfn(spikes(:,si));
         pk(cnt,:) = pk_;
      end
      si = si + 1;
   end
   
   % If at least 1 spike matched the template (hopefully 5 did)
   if cnt > 0
      mu  = mean(pk);
      sig = var(pk);
   else
      mu  = peakfn( fsp );
      sig = sqrt( max( fsp ) - min( fsp ) ) * ones(1,2);
   end
end













