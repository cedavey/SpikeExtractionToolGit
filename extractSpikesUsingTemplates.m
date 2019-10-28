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
   peakN     = round( mean( peakind ) );
   % peak function calculates the peak of the spike - e.g. using total diff
   % btwn min & max, or just max
   peakfn    = @(sp) sp(peakN) - min(sp); 
   peakfn    = @(sp) sp(peakN); 
   peakfn    = @(sp) [ min(sp) sp(peakN) ]; 
   % peakdiff function calculates diff btwn peaks, e.g. btwn total range or
   % sums diff btwn maxes & mins, & perhaps accounts for time btwn spikes
   % - weight positive peak more heavily than minimum peak
   peakdiff  = @(curr_peak, last_peaks) sum( abs((curr_peak - last_peaks) .* [1/2 2] ./ last_peaks), 2 );
   % require both +ve & -ve peaks to be within change threshold, and then
   % take spike with the min change
   peaktest  = @(curr_peak, last_peaks) abs( (curr_peak - last_peaks) ./ last_peaks );
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
   change_rate = params.ap_peak_change.value/100; % allowed rate of peak amplitude change (as percentage)
   newTemplates= params.allow_new_aps.value;      % allow new AP templates or not
   % Lambda for families
   lambda = params.ap_peak_change.value/100;
   kappa = params.kappa.value; % Times variance
      
   % Print progress
   str = sprintf( '\tSorting spikes...\n' );
   cprintf( 'Keywords', str );
   
   if ~strcmp('none',opts.debugOption)
      mup_vector = cell(nAP,1); mup_vector{1}{1} = [];
      mun_vector = cell(nAP,1); mun_vector{1}{1} = [];
      varp_vector = cell(nAP,1); varp_vector{1}{1} = [];
      varn_vector = cell(nAP,1); varn_vector{1}{1} = [];
   end
   
   nS = size(spikes, 2);
   prevProg = 1;
   for si=1:nS % for each possible spike (format: time x spike)
      % Print current percentage
      if 10*floor(si/nS*10) ~= prevProg
          fprintf( '\t%s is approx %d%% done\n', method, 10*floor(si/nS*10) );
          prevProg = 10*floor(si/nS*10);
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
               % In order to have a match threshold independent of amplitude,
               % scale by the size of the template's signal. Since cov is
               % voltage squared, scale by variance rather than std dev
               rho = cov_x_Yvec(curr_spike, APtemplates) ./ var(APtemplates); 
               
         otherwise
               str = sprintf('%s is an unknown match type', matchtype);
               displayErrorMsg(str);
               return;
      end
      
      % If spike is sufficiently similar in shape to a template, 
      % add it to one of the template families, or start a new family 
      % if the spike amplitude is not sufficiently similar.
      if any(rho > matchthresh)
         k = compare_templates(rho, curr_spike, APtemplates); % template with closest match
         
         % if no spike families yet, create one
         if isempty(APfamilies{k})
            APfamilies{k}{1}.spikes    = curr_spike;
            APfamilies{k}{1}.stimes    = curr_stime;
            APfamilies{k}{1}.meanspike = curr_spike;
            APfamilies{k}{1}.last_peak = peakfn(curr_spike); % size of last spike in family
            APfamilies{k}{1}.last_time = timefn(curr_stime); % time of last spike
            % Distribution of peaks (mean and var)
            APfamilies{k}{1}.mup = APfamilies{k}{1}.last_peak(2); % Positive peak mean
            APfamilies{k}{1}.mun = APfamilies{k}{1}.last_peak(1); % Negative peak mean
            APfamilies{k}{1}.varp = 0; % Positive peak variance
            APfamilies{k}{1}.varn = 0; % Negative peak variance
            
            if ~strcmp('none',opts.debugOption)
               APfamilies{k}{1}.mup_vector = [APfamilies{k}{1}.mup];
               APfamilies{k}{1}.mun_vector = [APfamilies{k}{1}.mun];
               APfamilies{k}{1}.varp_vector = [APfamilies{k}{1}.varp];
               APfamilies{k}{1}.varn_vector = [APfamilies{k}{1}.varn];
            end
            
         % if we have spike families see if we're at the right
         % scale, so get size of current peak, & time diff btwn
         % last spike and this spike, & check if the peak is
         % sufficiently the same to be considered part of the fam
         else
            curr_peak     = peakfn(curr_spike);
            curr_time     = timefn(curr_stime);
            last_peaks    = getStructFieldFromCell(APfamilies{k}, 'last_peak', 'array');
            last_time     = getStructFieldFromCell(APfamilies{k}, 'last_time', 'array');
            
            varn          = getStructFieldFromCell(APfamilies{k}, 'varn', 'array');
            varp          = getStructFieldFromCell(APfamilies{k}, 'varp', 'array');
            mun           = getStructFieldFromCell(APfamilies{k}, 'mun', 'array');
            mup           = getStructFieldFromCell(APfamilies{k}, 'mup', 'array');
            
            meanspike = getStructFieldFromCell(APfamilies{k}, 'meanspike', 'array');
                        
            % get rate of change from last peak (this is populated with the
            % mean spike now) to current peak, but then switch around
            % because I was finding that spikes were always shrinking 
            % over time, because if last_peaks is smaller (as in we're 
            % growing) when we divide by last peaks you're dividing by a 
            % small number, which makes the change rate larger, & less
            % likely to be within the allowed limit
%             change_rates1 = peakdiff( curr_peak, last_peaks );
%             change_test1  = peaktest( curr_peak, last_peaks );
%             change_rates2 = peakdiff( last_peaks, curr_peak );
%             change_test2  = peaktest( last_peaks, curr_peak );
            change_rates1 = peakdiff( curr_peak, [mun mup] );
            change_test1  = peaktest( curr_peak, [mun mup] );
            change_rates2 = peakdiff( [mun mup], curr_peak );
            change_test2  = peaktest( [mun mup], curr_peak );
%             mspk = zeros(size(meanspike,2), 2);
%             for ii = 1:size(meanspike,2), mspk(ii,:) = peakfn(meanspike(:,ii)); end
%             change_rates1 = peakdiff( curr_peak, mspk );
%             change_test1  = peaktest( curr_peak, mspk );
%             change_rates2 = peakdiff( mspk, curr_peak );
%             change_test2  = peaktest( mspk, curr_peak );
            
            % use the change_rates with the smallest value since it'll be
            % the min change rate that's selected below
            change_rates  = ternaryOp( min( change_rates1(:) ) < min( change_rates2(:) ), ...
                                       change_rates1, change_rates2 );
                                    
            % allow the spike to grow or shrink - merge so either direction
            % can pass the validity test
            change_test   = [ change_test1; change_test2 ];
            % test with change rate / 2 since it can go above or below the
            % mean spike for the family. 
            % * Wasn't really doing much with the /2, so I removed it.
%             valid1        = all( change_test1 < change_rate, 2 ); % all( change_test1 < change_rate/2, 2 ); 
%             valid2        = all( change_test2 < change_rate, 2 ); % all( change_test2 < change_rate/2, 2 ); 
%             valid         = valid1 | valid2;
            
            % Change the validation method to within range of a gaussian
            % distribution's mu and variance
            stdn_ = sqrt(varn);
            stdp_ = sqrt(varp);
            stdn_(varn == 0) = change_rate * abs(mun(varn == 0) / kappa);
            stdp_(varp == 0) = change_rate * abs(mup(varp == 0) / kappa);
            valid1 = zeros(size(mup));
            valid2 = zeros(size(mun));
            for j = 1:numel(mup)
               valid1(j) = (curr_peak(1) >= mun(j) - kappa * stdn_(j)) & (curr_peak(1) <= mun(j) + kappa * stdn_(j));
               valid2(j) = (curr_peak(2) >= mup(j) - kappa * stdp_(j)) & (curr_peak(2) <= mup(j) + kappa * stdp_(j));
            end
            valid = valid1 | valid2;
                        
            if any(valid)
               % Adding one spike to current axon
               [minrate, mini] = min( change_rates(valid) );
               % if current spike peak is within required rate of
               % change in peak amplitude for a spike family, add it
               % to the fam. Else start a new family.
               valid_ind = find( valid );
               i = valid_ind( mini );
               APfamilies{k}{i}.spikes(:,end+1) = curr_spike;
               APfamilies{k}{i}.stimes(:,end+1) = curr_stime;
               APfamilies{k}{i}.meanspike       = mean( APfamilies{k}{i}.spikes, 2 ); 

               % Distribution of peaks (mean and var)
               N_ = size(APfamilies{k}{i}.spikes,2);
               mup_prev = APfamilies{k}{i}.mup;
               mun_prev = APfamilies{k}{i}.mun;
               
%                mup_curr = (APfamilies{k}{1}.last_peak(2) + ((N_ - 1) * mup_prev * lambda))/N_; % Positive peak mean
               mup_curr = mup_prev + (APfamilies{k}{i}.last_peak(2) - mup_prev); % Positive peak mean
               mup_curr = mup_prev * lambda + mup_curr * (1 - lambda);
               APfamilies{k}{i}.mup = mup_curr;
               
%                mun_curr = (APfamilies{k}{1}.last_peak(1) + ((N_ - 1) * mun_prev * lambda))/N_; % Negative peak mean
               mun_curr = mun_prev + (APfamilies{k}{i}.last_peak(1) - mun_prev); % Negative peak mean 
               mun_curr = mun_prev * lambda + mun_curr * (1 - lambda);
               APfamilies{k}{i}.mun = mun_curr;
               
%                Snp_prev = (N_-1)*(APfamilies{k}{i}.varp^2); % Previous standard deviation
%                Snn_prev = (N_-1)*(APfamilies{k}{i}.varn^2); % Previous standard deviation
%                APfamilies{k}{i}.varp = abs( sqrt((Snp_prev + (APfamilies{k}{i}.last_peak(2) - mup_prev) * (APfamilies{k}{i}.last_peak(2) - mup_curr)) / N_) ); % Positive peak variance
%                APfamilies{k}{i}.varn = abs( sqrt((Snn_prev + (APfamilies{k}{i}.last_peak(1) - mun_prev) * (APfamilies{k}{i}.last_peak(1) - mun_curr)) / N_) ); % Negative peak variance
               varp_prev = APfamilies{k}{i}.varp;
               varn_prev = APfamilies{k}{i}.varn;
               varp_curr = varp_prev + ((APfamilies{k}{i}.last_peak(2) - mup_curr) .* (APfamilies{k}{i}.last_peak(2) - mup_prev) - varp_prev);
               varn_curr = varn_prev + ((APfamilies{k}{i}.last_peak(1) - mun_curr) .* (APfamilies{k}{i}.last_peak(1) - mun_prev) - varn_prev);
               APfamilies{k}{i}.varp = ternaryOp(varp_prev == 0, varp_curr, varp_prev * lambda  +  varp_curr * (1-lambda)); % Avoid making var small if it's only the second value
               APfamilies{k}{i}.varn = ternaryOp(varn_prev == 0, varn_curr, varn_prev * lambda  +  varn_curr * (1-lambda)); % Avoid making var small if it's only the second value
               
               if ~strcmp('none',opts.debugOption)
                  APfamilies{k}{i}.mup_vector = [ APfamilies{k}{i}.mup_vector mup_curr];
                  APfamilies{k}{i}.mun_vector = [ APfamilies{k}{i}.mun_vector mun_curr];
                  APfamilies{k}{i}.varp_vector = [ APfamilies{k}{i}.varp_vector APfamilies{k}{i}.varp];
                  APfamilies{k}{i}.varn_vector = [ APfamilies{k}{i}.varn_vector APfamilies{k}{i}.varn];
               end
               
               % getting peak to peak from family mean now, because it
               % wasn't working at all well using the last spike
               APfamilies{k}{i}.last_peak  = peakfn( curr_spike ); % size of last spike in family
%                APfamilies{k}{i}.last_peak   = peakfn( APfamilies{k}{i}.meanspike ); % size of last spike in family
               APfamilies{k}{i}.last_time   = timefn( curr_stime ); % time of last spike
            else
               % Creating new axon
               APfamilies{k}(end+1)         = cell(1);
               APfamilies{k}{end}.spikes    = curr_spike;
               APfamilies{k}{end}.stimes    = curr_stime;               
               APfamilies{k}{end}.meanspike = curr_spike; 
               APfamilies{k}{end}.last_peak = peakfn( curr_spike ); % size of last spike in family
               APfamilies{k}{end}.last_time = timefn( curr_stime ); % time of last spike
                              
               APfamilies{k}{end}.mup = APfamilies{k}{end}.last_peak(2); % Positive peak mean
               APfamilies{k}{end}.mun = APfamilies{k}{end}.last_peak(1); % Negative peak mean
               APfamilies{k}{end}.varp = 0; % Positive peak variance
               APfamilies{k}{end}.varn = 0; % Negative peak variance
               
               if ~strcmp('none',opts.debugOption)
                  APfamilies{k}{end}.mup_vector = [ mup_curr];
                  APfamilies{k}{end}.mun_vector = [ mun_curr];
                  APfamilies{k}{end}.varp_vector = [APfamilies{k}{end}.varp];
                  APfamilies{k}{end}.varn_vector = [APfamilies{k}{end}.varn];
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
            % Forget factor. lambda defines how much weight do older spikes
            % have
            lambda = 0.95;
            APtemplates(:,k) = APtemplates(:,k) * lambda + (1 - lambda) * curr_spike;
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
end
















