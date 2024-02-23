% [APspikes,APtimes,AP_alignment_inds] = extractSpikesUsingTemplates( APtemplates, APnumsamples, AP_alignment_inds, tseries, method, params, normAPs ,opts)
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
%  params      - parameters for extrac14/9/2023ting spikes, set using SET gui
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
%
% 
function [APspikes, APtimes,AP_alignment_inds] = extractSpikesUsingTemplates( APtemplates, APnumsamples, AP_alignment_inds, tseries, method, params, normAPs ,opts)
  
% preallocating and doing settings
   if exist('opts', 'var') && isfield(opts,'loadingWindowOn'), progress_window = opts.loadingWindowOn; 
   else, progress_window = false; end
   if strcmpi( opts.auto_params, 'true' )||(isnumeric(opts.auto_params) && opts.auto_params==1)
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
%       peakfn    = @(sp) max(sp) - min(sp); 
%       peakfn    = @(sp) max(sp); 
      peakfn    = @(sp) [ min(sp) max(sp) ]; 
   else
      uniqueAPLength = false; 
      [nT, nAP] = size(APtemplates); % num samples in each template, & num templates
      peakind   = getMaxInd( APtemplates, 1 );
      peakN     = round( median( peakind ) ); % mean index of template peaks
      % peak function calculates the peak of the spike - e.g. using total diff
      % btwn min & max, or just max
%       peakfn    = @(sp) sp(peakN) - min(sp); 
%       peakfn    = @(sp) sp(peakN); 
      peakfn    = @(sp) [ min(sp) sp(peakN) ]; 
   end

   APspikes  = [];  APtimes = [];

   % peakdiff function calculates diff btwn peaks, e.g. btwn total range or
   % sums diff btwn maxes & mins, & perhaps accounts for time btwn spikes
   % and expects inputs of the form [peak_neg  peak_pos]
   % - weight positive peak more heavily than minimum peak
	% peakdiff  = @(curr_peak, last_peaks) sum( abs((curr_peak - last_peaks) .* [1/2 2] ./ last_peaks), 2 );
	% peakdiff  = @(x, y) abs(((x(:,2)-x(:,1))-(y(:,2)-y(:,1)))./(y(:,2)-y(:,1))) ;

   % divide by 2 to avg diff across pos & neg peaks
   peakdiff  = @(curr_peak, last_peaks) sum( abs((curr_peak - last_peaks) .* [1 1] ./ last_peaks), 2 ) / 2;
   timefn    = @(sp, st) st( getMaxInd(sp) );
% FIND ALL SPIKES
   switch lower(method)
      case 'matched filter'
         if uniqueAPLength
            [spikes, stimes, spikesFull, stimesFull,alignment_inds] = getSpikesByThresholding( tseries, params );
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
         str = sprintf('Error getting spikes by thresholding. Method %s is not implemented. Check ''extractSpikesUsingTemplates.m''.', method);
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
   % if kappa are 0, then allow any amount of change
   if kappa_neg==0, kappa_neg=Inf; end
   if kappa_pos==0, kappa_pos=Inf; end
      
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
            % use lagged corr in case the spike peak causes delay/advance
            matchfn = @( sp, temp ) max(xcorr(temp, sp, 'normalized', 5)); % consider max 5 lags

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
   % backup noise if getting a white block of samples fails (taken from 
   % the start of the signal acquisition so doesn't change with expmt changes)
   noise_backup  = getNoise(stimes, 1, tseries); 
   % keep a list of new axons so we can see if it was generated from a
   % stray spike, and actually becomes very similar to another axons quickly
   new_axon_check= 40;
   num_new_axons = 0; 
   is_new_axon   = false;
   progress_str  = []; % allows us to write over prev progress rather than a new line every time it's report 
   try
      for si=1:nS % for each possible spike (format: time x spike)
         
         
         [prevProg,progress_str] = print_progress_percentage(si,nS,prevProg,method,progress_str);
        

         % Update Progress window
         if progress_window
            try
               waitbar_handles = waitbar(si/nS, waitbar_handles,...
                  [method ' is ', num2str(round( si/nS*100 )), '% done']);
            catch
               delete(waitbar_handles);
               str = ['The process has been manually stopped on time: ', num2str(stimes{si}(1)),' seconds.\n'];
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
         rho = spikeCorrWithTemplates( curr_spike, APtemplates,alignment_inds(si),AP_alignment_inds, matchfn );

         % If spike is sufficiently similar in shape to a template, 
         % add it to one of the template families, or start a new family 
         % if the spike amplitude is not sufficiently similar.
         % - gotta account for cov so rho has to be within band around 1
         if any( rho > matchthresh )
            k = compare_templates(rho, curr_spike, APtemplates, APfamilies,alignment_inds(si), matchthresh); % template with closest match

            % if no spike families yet, create one
            if isempty(APfamilies{k})
                
               %ADD A NEW AXON TO THE FAMILY
               [APfamilies{k}{1},is_new_axon(k,1),num_new_axons(k,1)] = add_new_axon_to_family(stimes,si,curr_spike,alignment_inds(si),tseries,noise_backup,peakfn,timefn,opts);
            else
                
                
               %There are axons within the family, so let's see which one
               %fits
               [axon_ind] = match_spike_to_axon_using_peaks(APfamilies{k},curr_spike);
               
               if axon_ind ~= 0
                  %WE FOUND A GOOD MATCH, ADD SPIKE TO AXON
                  was_new = is_new_axon(k,axon_ind);
                  [APfamilies{k}{axon_ind},num_new_axons(k,axon_ind),is_new_axon(k,axon_ind)] = update_axon_with_spike(APfamilies{k}{axon_ind},curr_spike,curr_stime,alignment_inds(si),uniqueAPLength,is_new_axon(k,axon_ind),num_new_axons(k,axon_ind));
                 % if we have enough axons for a decent mean, check if
                 % it's very similar to another family, & merge
                 if was_new == true && is_new_axon(k,ai)==false
                    [APfamilies, is_new_axon, num_new_axons] = ...
                             mergeSimilarAxonFamilies( APfamilies, [k ai], is_new_axon, num_new_axons, ...
                                                       matchthresh, kappa_pos, kappa_neg );
                 end
                 
               else
                  % There was no good match - add a new axon to the family
                  %ADD A NEW AXON TO THE FAMILY
                  [APfamilies{k}{end+1},is_new_axon(k,end+1),num_new_axons(k,end+1)] = add_new_axon_to_family(stimes,si,curr_spike,alignment_inds(si),tseries,noise_backup,peakfn,timefn,opts);
               end
            end

%             % don't update AP templates as they're typically already
%             % comprised of several hundred samples, so adding a few more
%             % won't make much difference
%             % update template by incorporating current spike into
%             % mean - we don't care about amplitue
%             if normAPs
%                curr_spike  = curr_spike / std(curr_spike); 
%             end
%             % Equal importance to all past spikes
%             % APtemplates(:,k) = (APtemplates(:,k)*APnumsamples(k) + ...
%             % curr_spike) / (APnumsamples(k)+1); % uncomment this line and
%             % comment the next one to go back
%             % Forgetting factor - lambda defines how much weight do older spikes
%             % have
%             if uniqueAPLength
%                curr_spike = alignSpikeToTemplate( APtemplates{k}, curr_spike ); 
%                APtemplates{k} = APtemplates{k} * lambda  +  curr_spike * (1 - lambda);
%             else
%                APtemplates(:,k) = APtemplates(:,k) * lambda  +  curr_spike * (1 - lambda);
%             end
%             APnumsamples(k)  = APnumsamples(k) + 1;

         % Spike is not sufficiently similar to one of the families
         % so start a new AP template if user has allowed new ones
         else % rho does not exceed threshold
            if newTemplates
               % make a new AP family within the templates
               if uniqueAPLength 
                  APtemplates{end+1} = curr_spike;
               else
                  APtemplates(:, end+1) = curr_spike;
               end
               APnumsamples(end+1)   = 1;
               % Creating new axon family
               % Get noise samples prior to the current spike to initialize
               % the variance from the noise.
               [APfamilies{end+1}{1},is_new_axon(end+1,1),num_new_axons(end+1,1)] = add_new_axon_to_family(stimes,si,curr_spike,alignment_inds(si),tseries,noise_backup,peakfn,timefn,opts);
            else
               bin_spikes(si) = 1; 
            end
         end % end if rho > match_thresh
      end % end for each spike
   catch ME
      catchstr = getCatchMEstring( ME, 'main: ' );
%        runtimeErrorHandler(ME);
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
      for ai=1:length( APfamilies )
         if ~isempty(APfamilies{ai})
             nf = cellfun( @(fam) size( fam.spikes,2 ), APfamilies{ai} );
             toosmall = nf <= params.min_spiking_threshold.value; 
             APfamilies{ai}( toosmall ) = [];
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
         nP = 4; % number of features to plot (+/- mean & +/- var)
         if nAP<4
            nc = nAP; nr = nP;
         else
            % too many templates --> just plot in pairs
            nc = ceil( sqrt( nAP*nP ) ); nc = ternaryOp( round(nc/2)~=nc/2, nc+1, nc );
            nr = ceil( nAP*nP / nc );
         end
         Nfam = length(APfamilies);
         for ti = 1:Nfam % for each template
            mu_pos = -Inf; var_pos = -Inf;
            mu_neg =  Inf; var_neg = -Inf;
            naxons = length( APfamilies{ti} );
            cols   = getColourMatrix( naxons );

            lw = 2;
            for ai = 1:naxons % for each family within the template
               % if unique length spikes then NaNs exist where spike
               % doesn't have a value cuz it's shorter than others -
               % therefore construct time vector based position of peak
               % spike value within each spike
               sp        = APfamilies{ti}{ai}.spikes(:,1);
               [~,pind]  = nanmax(sp);
               tvec      = APfamilies{ti}{ai}.stimes(pind,:);
               peakvec   = APfamilies{ti}{ai}.spikes(pind,:); % aligned at peak
               troughvec = nanmin( APfamilies{ti}{ai}.spikes, [], 1 );

               subplot(nr,nc,0*nAP+ti); hold on;
                  plot( tvec, APfamilies{ti}{ai}.mup_vector, 'color', cols(ai,:) );
                  legend('\mu_{pos}');

               subplot(nr,nc,1*nAP+ti); hold on;
                  plot( tvec, APfamilies{ti}{ai}.mun_vector, 'color', cols(ai,:) );
                  legend('\mu_{neg}');

               subplot(nr,nc,2*nAP+ti); hold on;
                  plot( tvec, +APfamilies{ti}{ai}.varp_vector, '-',  'color', cols(ai,:), 'linewidth', lw );
                  plot( tvec, -APfamilies{ti}{ai}.varn_vector, '--', 'color', cols(ai,:), 'linewidth', lw );
                  var_pos = ternaryOp( var_pos < max( APfamilies{ti}{ai}.varp_vector ), max( APfamilies{ti}{ai}.varp_vector ), var_pos );
                  var_neg = ternaryOp( var_neg < max( APfamilies{ti}{ai}.varn_vector ), max( APfamilies{ti}{ai}.varn_vector ), var_neg );
                  legend('\sigma^2_{pos}', '\sigma^2_{neg}');

               subplot(nr,nc,3*nAP+ti); hold on;
                  plot( tvec, peakvec,   '-', 'color', cols(ai,:), 'linewidth', lw );
                  plot( tvec, troughvec, '-', 'color', cols(ai,:), 'linewidth', lw );
                  legend('spike_{max}', 'spike_{min}');

            end
         end

         % now plot all family & axon peaks on the same plot
         figure; 
         Faxons = cellfun(@length, APfamilies); % number of axons within each family
         Naxons = sum(Faxons); % total number of distinct axons
         cols   = getColourMatrix( Naxons );

         lw = 2; pp = 1; % number of peak plots so far
         legstr = cell(0);
         for ti = 1:length(APfamilies) % for each template
            naxons = length( APfamilies{ti} ); % axons within this family

            for ai = 1:naxons % for each family within the template
               sp        = APfamilies{ti}{ai}.spikes(:,1);
               [~,pind]  = nanmax(sp);
               tvec      = APfamilies{ti}{ai}.stimes(pind,:);
               peakvec   = APfamilies{ti}{ai}.spikes(pind,:); % aligned at peak
               troughvec = nanmin( APfamilies{ti}{ai}.spikes, [], 1 );
               legstr{end+1} = sprintf('fam %d, axons %d', ti, ai);

               subplot(311); hold on;
                  plot( tvec, peakvec,   '.', 'color', cols(pp,:) );
                  title('Axon peaks');

               subplot(312); hold on;
                  plot( tvec, troughvec, '.', 'color', cols(pp,:) );
                  title('Axon troughs');

               subplot(313); hold on;
                  plot( tvec, peakvec,   '.', 'color', cols(pp,:) );
                  plot( tvec, troughvec, '.', 'color', cols(pp,:) );
                  title('Axon peaks and troughs');

               pp = pp+1;
            end
         end
         a1 = subplot(311);
            yl = get(a1,'ylim'); set(a1,'ylim', [0 yl(2)]);
         a2 = subplot(312);
            yl = get(a2,'ylim'); set(a2,'ylim', [yl(1) 0]);
         if ~isempty(legstr), legend(legstr); end
      end

    
      % Each valid spike has been allocated to an AP template, & to an AP
      % familiy within that. Now for each AP family we need to construct a
      % timeseries
      [APspikes,APtimes] = construct_APfamily_timeseries(APfamilies,nAP,params,tseries);
      
      
   catch ME
      str = getCatchMEstring( ME, 'end: ' );
%       runtimeErrorHandler(ME);
   end
end
% Print current percentage - old progress text is deleted from command
% window and new added 
function [prevProg,progress_str] = print_progress_percentage(si,nS,prevProg,method,progress_str);
         
         if 10*floor(si/nS*10) ~= prevProg
            str = sprintf( '\n\t%s is approx %d%% done\n', method, 10*floor(si/nS*10) );
            refreshdisp(str, progress_str);
            progress_str = str;
            prevProg = 10*floor(si/nS*10);
         end
end

function [mu, sig] = initialize_mu(curr_spike, noise_samples)
  
   mu  = [min(curr_spike) max(curr_spike)];    
   sig = var(noise_samples);   % variance of noise
   sig = [sig sig];            % variance for peak & trough of spike
    
    % Before I was extracting the first N spikes that were correlated with
    % the input spike, but this doesn't ensure that the spikes are the
    % right scale, so it was fucking up the mean & variance estimates.
    % Instead just use the first spike, but allow distance from the mean to
    % be larger than the user requested number of std dev's from the mean,
    % until there is a sufficient number of spikes to have converged to the mean 
%    maxlag = 2; % number of lags in cross-corr
%    % Initialize the mean peak for a new family, triggered by spike at index 
%    % si not having a family of similar size. Goes through the first N
%    % spikes that match the template and chooses that value as the mean.
%    if iscell( spikes ), uniqueAPLength = true; else, uniqueAPLength = false; end
%    if uniqueAPLength
%       nSp = length( spikes );
%       fsp = spikes{si}; % spike we're making the family from
%    else
%       nSp = size( spikes, 2 );
%       fsp = spikes(:,si);
%    end
%    fam_peaks = [ min(fsp) max(fsp) ];
%    si  = si + 1; % start matching with spikes after current spike
%    allowed_diff = 0.5; % Allowed difference between spikes to initialize mu
%    % if we've run outta spikes in the spike train, return
%    if si >= nSp
%       mu = peakfn( fsp ); sig = 0; 
%       return; 
%    end
%    
%    cnt = 0;
%    pk = [];
%    while cnt < 10 && si < nSp
%       if uniqueAPLength
%          sp = spikes{si}; % spike we're making the family from
%       else
%          sp = spikes(:,si);
%       end
%       sp_peaks = [ min(sp) max(sp) ];
%       % align the spikes at their peaks 
%       [ sp, fsp ] = alignDiffLengthSpikes( sp, fsp );
%       c   = xcorr( sp, fsp, maxlag, 'coeff' );
%       rho = max( abs(c) );
%       d1  = diffpeakfn( sp_peaks, fam_peaks ); % diff in peaks btwn family spike & sp
%       d2  = diffpeakfn( fam_peaks, sp_peaks ); % diff in peaks in reverse order
%       d   = min( d1, d2 ); % get largest diff btwn pos & neg peak changes
% %       d   = max( d1, d2 ); % get largest diff btwn pos & neg peak changes
% 
%       % If spike is closeley correlated to the template
%       if rho > matchthresh && d < max(fsp)*0.3 % TO DO: DON'T HARD-CODE ALLOWABLE DIFFERENCE !! 
%          cnt     = cnt + 1;
%          if uniqueAPLength
%             pk_  = peakfn( spikes{si} );
%          else
%             pk_  = peakfn(spikes(:,si));
%          end
%          pk(cnt,:) = pk_;
%       end
%       si = si + 1;
%    end
%    
%    % If at least 2 spikes matched the template (hopefully 5+ did!)
%    if cnt > 2
%        mu  = mean(pk, 1);
%        
%        if noise_samples == 0
%            % We are setting sig from the noise now. This line only works if no
%            % noise samples could be found
%              sig = var(pk, [], 1 ); 
%        else
%            sig = var(noise_samples) * ones(1,2);
%        end
%    else
%        if uniqueAPLength
%            fsp = spikes{:,si};
%        else
%            fsp = spikes(:,si);
%        end
%            
%        mu  = peakfn( fsp );
%        if noise_samples == 0
%            % We are setting sig from the noise now. This line only works if no
%            % noise samples could be found
%            sig = sqrt( max( fsp ) - min( fsp ) ) * ones(1,2);
%        else
%            sig = var(noise_samples) * ones(1,2);
%        end
%    end
end

% Extract some noise samples from a few samples before the current spike's
% onset.
function noise_samples = getNoise(stimes, si, tseries, noise_backup)

    Nsamples  = 400;
    Nmin      = 20;
    havenoise = false; 

    % try getting 100 samples of white noise - if not white reduce size
    % until find a contiguous block of signal that is white noise
    current_sample = ceil(stimes{si}(1) / tseries.dt);
    idx            = [max(current_sample - Nsamples, 1):current_sample];
    noise_samples  = tseries.data(idx);
    while ~havenoise && idx(1)>1 && Nsamples>Nmin
        % test current noise sample to see if it's white
        if iswhite(noise_samples)
            havenoise = true; % successfully have contiguous block of white noise
            return;
            
        % reduce size of contiguous block to see if we can get white noise
        else
           % trim from both ends just in case we have a spike right at the
           % beginning or end, so end up with no white noise
            Ntrim          = round(Nsamples*1/3);     % trim 1/3, keep 2/3
            Ntrim          = ceil(Ntrim/2);           % trim from both ends
            Nsamples       = round(Nsamples * 2/3);
            idx            = [max(current_sample - Nsamples - Ntrim, 1):(current_sample - Ntrim)];
            noise_samples  = tseries.data(idx);
            current_sample = idx(end);                % shift current sample by Ntrim
        end
    end
    
    % If we're here we were unsuccessful in finding a contiguous block of
    % white noise samples that has length at least Nmin. Try setting length
    % and moving back through the timeseries, rather than staying where we
    % are & reducing the length of the timeseries
    if Nsamples < Nmin
        Nsamples       = 100; % 50 samples to check for noise
        current_sample = ceil(stimes{si}(1) / tseries.dt);
        idx            = [max(current_sample - Nsamples, 1):current_sample];
        noise_samples  = tseries.data(idx);
        while idx(1) > 1
            if iswhite(noise_samples)
                % If white, return the noise
                return
            else
                % If not white, take previous samples
                idx = idx - Nsamples;
                noise_samples = tseries.data(idx(idx>0));
            end
        end
        
        % If everything failed, re-start with index 'si' and start
        % reducing the sample size one by one until size = Nmin, then start
        % moving forward. (Artemio 29/Nov/2022)
        Nsamples = 100;
        current_sample = ceil(stimes{si}(1) / tseries.dt);
        while Nsamples > Nmin
            idx = [current_sample : min(current_sample + Nsamples, length(tseries.data))];
            noise_samples  = tseries.data(idx);
            if ~iswhite(tseries.data(idx))
                Nsamples = Nsamples-1;
            else
                return
            end
        end
        % If still fails, rather than reducing the size of the contiguous 
        % block, we try shifting forward to find a block of white noise
        current_sample = ceil(stimes{si}(1) / tseries.dt);
        Nsamples = 30;
        idx = [current_sample : min(current_sample + Nsamples, length(tseries.data))];
        while idx(end) <  length(tseries.data)
            noise_samples  = tseries.data(idx);
            if ~iswhite(tseries.data(idx))
                current_sample = current_sample + 1;
                idx = [current_sample : min(current_sample + Nsamples, length(tseries.data))];
            else
                return
            end
        end
        
        
        % If it reaches this point, it didn't find any noisy period
        if exist('noise_backup', 'var')
            noise_samples = noise_backup;
        else
            % If no approach found a long enough period of white noise and
            % there was no 'noise_backup' input to this function, return 
            % the first Nmin samples of the recording. It is rarely
            % expected for the code to enter this segment.
            noise_samples = tseries.data(1:Nmin);
        end
    end
end

% Spawn a new axon family with all the necessary struct fields
function family = newAxonFamily( curr_spike, peakfn, curr_stime, timefn, fam_mu, fam_var, alignment_ind, opts )
   family.spikes    = curr_spike;
   family.stimes    = curr_stime;
   family.meanspike = curr_spike;
   family.last_peak = peakfn( curr_spike ); % size of last spike in family
   family.last_time = timefn( curr_spike, curr_stime ); % time of last spike
   family.mun       = fam_mu(1); 
   family.mup       = fam_mu(2); 
   family.varn      = fam_var(1);
   family.varp      = fam_var(2);
   family.Nspikes   = 1;
   family.alignment_ind = alignment_ind;
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

%Assemble categorised spikes into a timeseries
function [APspikes,APtimes] = construct_APfamily_timeseries(APfamilies,nAP,params,tseries)
    [APspikes,APtimes] = deal(cell(1,nAP));
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
                  sind(isnan(sind)) = [];
%                   peakind   = max( fam{ff}.spikes, [], 1);
                  fam_tseries(sind,ff) = fam{ff}.spikes(valid);
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
end

%initialise a new axon within a family
function [APfam_axon,is_new_axon,num_new_axons] = add_new_axon_to_family(stimes,si,curr_spike,alignment_inds,tseries,noise_backup,peakfn,timefn,opts)
       % Distribution of peaks (mean and var)
       % Get noise samples prior to the current spike to initialize
       % the variance from the noise.
       noise_samples     = getNoise(stimes, si, tseries, noise_backup);
       % Initialize mean and var - can't use template since it
       % gives shape, but not scale
       [fam_mu, fam_var] = initialize_mu( curr_spike, noise_samples);
       APfam_axon = newAxonFamily( curr_spike, peakfn, curr_stime, timefn, fam_mu, fam_var, alignment_inds(si), opts );
       is_new_axon  = true; 
       num_new_axons= 1; % record spike for new axon
        % if we have spike families see if we're at the right
            % scale, so get size of current peak, & time diff btwn
            % last spike and this spike, & check if the peak is
            % sufficiently the same to be considered part of the fam
end


%If thare are any valid matches between the spike and the axons within
%APfam, the best one is selected (index of the best is output), else 0 is
%output.
function [axon_ind] = match_spike_to_axon_using_peaks(APfam,curr_spike)
      curr_peak = peakfn(curr_spike);           
               varn      = getStructFieldFromCell( APfam, 'varn', 'array' );
               varp      = getStructFieldFromCell( APfam, 'varp', 'array' );
               mun       = getStructFieldFromCell( APfam, 'mun', 'array' );
               mup       = getStructFieldFromCell( APfam, 'mup', 'array' );

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

               % Because APtemplate gives shape but not scale, the first
               % spike is used to initialise the axon's mean, which is a
               % noisy estimate of the mean. Therefore, allow the
               % distance from the mean to be larger at the beginning,
               % until we've got 20-ish spikes to get a decent estimate
%                Nconverge     = 20; % assume that mean est has converged by now
%                scaleStd      = 10; % extra percent of std devs we'll allow at the beginning
%                Naxons        = getStructFieldFromCell( APfamilies{k}, 'Nspikes', 'array' );
%                % add scaleStd percent to number of std devs from mean (kappa)
%                extra         = max( (Nconverge - Naxons) / Nconverge * scaleStd/100 + 1, 1 ); % e.g., additional 10%
               extra  = ones(numel(mup),1); % initial variance estimate is good now, so don't need the extra

               stdn_  = sqrt(varn);
               stdp_  = sqrt(varp);
               valid1 = zeros(size(mup));
               valid2 = zeros(size(mun));
               try
                  % calculate the probability of a spike belonging to a
                  % particular neuron using Bayesian inference.
                  % Prior probability: current firing rate rather than
                  % number of spikes to allow for axons to start and stop
                  % firing
                  for ti = 1:numel(mup)
                     trough_lowerbound = ( mun(ti) - (kappa_neg*extra(ti)) * stdn_(ti) );
                     trough_lowerbound = ternaryOp( trough_lowerbound<=0, trough_lowerbound, 0);

                     trough_upperbound = ( mun(ti) + (kappa_neg*extra(ti)) * stdn_(ti) );
                     trough_upperbound = ternaryOp( trough_upperbound<=0, trough_upperbound, 0);

                     peak_lowerbound   = ( mup(ti) - (kappa_pos*extra(ti)) * stdp_(ti) );
                     peak_lowerbound   = ternaryOp( peak_lowerbound>=0, peak_lowerbound, 0);
                     
                     peak_upperbound   = ( mup(ti) + (kappa_pos*extra(ti)) * stdp_(ti) );
                     peak_upperbound   = ternaryOp( peak_upperbound>=0, peak_upperbound, 0);
                     if trough_lowerbound > 0 || trough_upperbound > 0
                        str = sprintf('trough = %.2f, lower = %.2f, upper = %.2f\n', ...
                                      curr_peak(1), trough_lowerbound, trough_upperbound );
                        cprintf('Error*', str );
                     end
                     if peak_lowerbound < 0 || peak_upperbound < 0
                        str = sprintf('peak = %.2f, lower = %.2f, upper = %.2f\n', ...
                                      curr_peak(2), peak_lowerbound, peak_upperbound );
                        cprintf('Error*', str );
                     end

                     valid1(ti) = ( curr_peak(1) >= trough_lowerbound ) ...
                               && ( curr_peak(1) <= trough_upperbound );
                     valid2(ti) = ( curr_peak(2) >= peak_lowerbound ) ...
                               && ( curr_peak(2) <= peak_upperbound ); %
                  end
               catch ME
                  disp('err');
               end
               valid = valid1 & valid2;

               if any(valid)
                  % Adding one spike to axon family of best fit
                  [~, mini] = min( change_rates(valid) );
                  vind = find( valid );
                  axon_ind   = vind( mini );
               else
                  axon_ind = 0;
               end

end

%Now we've decided that this spike matches the axon, we need to the axon
%with the spike
function [Axon,num_new_axons,is_new_axon] = update_axon_with_spike(Axon,curr_spike,curr_stime,curr_alignment_ind,uniqueAPLength,is_new_axon,num_new_axons)
                  
               
      if uniqueAPLength

         [ Axon.meanspike, Axon.spikes, curr_spike, curr_stime, Axon.stimes,Axon.alignment_ind,~] = ...
               addSpikeToTemplateMaxCorr( Axon.meanspike, curr_spike,Axon.alignment_ind,curr_alignment_ind, Axon.spikes, curr_stime, Axon.stimes, lambda );
         Nspikes               = Axon.Nspikes + 1;
         Axon.Nspikes           = Nspikes;
      else
         Axon.Nspikes           = Axon.Nspikes + 1;
         Axon.spikes(:,Axon.Nspikes) = curr_spike;
         Axon.stimes(:,Axon.Nspikes) = curr_stime;
         % efficient calculation of meanspike for large families
         Axon.meanspike         = recursiveUpdate( Axon.meanspike , [], curr_spike,lambda );

      end
       
    try
      % Distribution of peaks (mean and var)
      [mup_curr, varp_curr] = recursiveUpdate( Axon.mup, Axon.varp, curr_peak(2), lambda );
      [mun_curr, varn_curr] = recursiveUpdate( Axon.mun, Axon.varn, curr_peak(1), lambda );
      Axon.mup       = mup_curr;
      Axon.mun       = mun_curr;
      Axon.varp      = varp_curr;
      Axon.varn      = varn_curr;             
      Axon.last_peak = peakfn( curr_spike ); % size of last spike in family
      Axon.last_time = timefn( curr_spike, curr_stime ); % time of last spike             
      if ~strcmp( 'none', opts.debugOption )
         Axon   = updateDebugInfo( Axon, mup_curr, mun_curr, varp_curr, varn_curr );
      end

      % check if this axon family is new, update its stats, and
      % merge with existing family if it turns out they're
      % heaps similar
                  

     catch ME
        disp('something is wrong');
     end
end
