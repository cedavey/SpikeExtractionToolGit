% generateSpikeStatistics(tseries, method, method_params)
% Generate statistics on rate timeseries for each AP family. Currently
% the only statistic is autocorrelation. 
function generateSpikeStatistics(tseries, method, method_params)
   maxax = 12; % max axons within each family to display
   nap   = length( tseries.data );
   nfam  = cellfun( @( d ) size(d,2), tseries.data );
   dt    = tseries.dt;
   
   switch lower(method)
      %% spike raster
      case 'raster'
         fopts = {'fontweight', 'bold', 'fontsize', 16};
         fh = figure; set(fh, 'Visible', 'off'); % figure for histogram of amplitudes
         
         % for each AP templates, plot the autocorrelation for each family 
         cols  = getColourMatrix(nap); % get each template a unique colour
         index = 0;
         for ii=1:nap
            [nT, nfam] = size( tseries.data{ii} ); 

            for jj=1:nfam
               stimes = round( tseries.APstimes{ii}{jj} );
               index  = index + 1;
               figure( fh );
               hold on;
               plot(stimes, index * ones(size(stimes)), '.', 'color', cols(ii,:) );
            end
         end
         set(fh, 'Visible', 'on'); % figure for histogram of amplitudes
%          set(gca, 'xcolor', [1 1 1]);
         set(gca, 'ycolor', [1 1 1]);
         xlabel('Time (seconds)', fopts{:});
         xlim( [tseries.time(1) tseries.time(end)] )
         ylim([1 index]);
      
      case 'psth'
         fopts = {'fontweight', 'bold', 'fontsize', 12};
         fh = figure; set(fh, 'Visible', 'off'); % figure for histogram of amplitudes
         numsecs  = method_params.time_post_spike.value;
         num_bins = method_params.num_bins.value;
         numsamp  = round( numsecs / dt );
         dbin     = numsecs / num_bins;
         timebins = dbin:dbin:(numsecs+dbin); % add list bin to catch everything above

         % for each AP templates, plot the autocorrelation for each family 
         cols = getColourMatrix(nap); % get each template a unique colour
         nc   = min( max( nfam ), 6 ); % max of 2 rows per AP family, max 6 cols per row
         nr   = sum(cellfun(@(d) min([2 ceil(size(d,2)/nc)]), tseries.data));         

         for ii=1:nap
            if nfam(ii) > 2*nc % max of 2 rows per family
               str  = sprintf('Only displaying %d from template family of %d', ...
                               maxax, nfam(ii));
               displayErrorMsg(str);
               nfam(ii) = 2*nc;
            end
            % start new APs axes at beginning of new row
            fi = ternaryOp( ii==1, 1, sum(ceil(nfam(1:ii-1)/nc))*nc + 1 );

            for jj=1:nfam(ii)
               % for each spike, count the number of spikes in each time
               % bins after it, up to numsecs after the spike
               stimes    = tseries.APstimes{ii}{jj}; % spike times list
               timecount = zeros(1, num_bins+1 ); % count for each time bin
               while length(stimes)>1
                  sptime = stimes(1);
                  stimes(1) = [];
                  postsp    = stimes( (stimes-sptime) <= numsecs ) - sptime;
                  count     = hist( postsp, timebins );
                  timecount = timecount + count; 
               end
               timecount = timecount(1:end-1);
               timecount = timecount / sum(timecount); % convert to prob
               figure(fh);
               subplot(nr, nc, fi), fi=fi+1;
%                   plot(timebins(1:end-1), timecount, '.', 'color', cols(ii,:) );
                  bar(timebins(1:end-1), timecount );
                  set(get(gca,'children'),'facecolor',cols(ii,:));
                  box off;
                  xlim([timebins(1)-(timebins(2)-timebins(1)) timebins(end)]);
                  xlabel('Time (seconds)', fopts{:});
                  xlabel('P(spike)', fopts{:});
            end
         end
         set(fh, 'Visible', 'on'); % figure for histogram of amplitudes
%          set(gca, 'xcolor', [1 1 1]);
         
      %% ISI
      case 'interspike interval'
         % for each AP templates, plot the ISI for each family          
         ft = figure; set(ft, 'Visible', 'off');
         num_bins = method_params.num_bins.value;
         max_isi  = method_params.max_isi.value;

         nc = 6; % max of 2 rows per AP family
         maxfam = max( cellfun(@(d) size(d,2), tseries.data) );
         nc = min( maxfam, nc );
         nr = sum(cellfun(@(d) min([2 ceil(size(d,2)/nc)]), tseries.data));         

         nfam = zeros(nap,1); % keep track of number of axes displayed
         for ii=1:nap
            [~, nfam(ii)] = size(tseries.data{ii}); 
            if nfam(ii) > 2*nc % max of 2 rows per family
               str  = sprintf('Only displaying %d from template family of %d', ...
                               maxax, nfam(ii));
               displayErrorMsg(str);
               nfam(ii) = 2*nc;
            end
            % start new APs axes at beginning of new row
            if ii==1
               fi = 1;
            else
               fi = sum(ceil(nfam(1:ii-1)/nc))*nc + 1;
            end
            
            max_count = 0;
            for jj=1:nfam(ii)
               stimes = tseries.APstimes{ii}{jj};
               isi    = diff(stimes);
               isi(isi>max_isi) = [];
               [count, bins] = hist(isi, num_bins);
               % figure(fa); subplot(nr, nc, fi); autocorr(isi,'NumLags',size(isi,2)-1);
               count  = count / sum(count);
               if max(count) > max_count, max_count = max(count); end
               str    = sprintf('AP %d, axon %d: %d spikes', ii, jj, length(stimes));
               figure(ft);
               subplot(nr, nc, fi), fi=fi+1;
                  bar(bins, count);
                  title(str);
                  xlabel('Time (s)');
                  ylabel('Probability');               
            end
            
            % Adjust y axis
            % start new APs axes at beginning of new row
            if ii==1
               fi_ = 1;
            else
               fi_ = sum(ceil(nfam(1:ii-1)/nc))*nc + 1;
            end
            for jj=1:nfam(ii)
               figure(ft);
               subplot(nr, nc, fi_), fi_=fi_+1;
               ylim([0 max_count]);
            end
         end
         set(ft, 'Visible', 'on');
         
      %% amplitude change over time
      case 'amplitude change'
         % for each AP templates, plot the autocorrelation for each family 
         ft = figure; set(ft, 'Visible', 'off'); % figure for amp change over time
         fh = figure; set(fh, 'Visible', 'off'); % figure for histogram of amplitudes
         mv =  Inf;
         Mv = -Inf;

         allpeaks = []; alltimes = [];
         for ii=1:nap
            [nT, nfam] = size(tseries.data{ii}); 
            if nfam > maxax
               str  = sprintf('Only displaying %d from template family of %d', ...
                               maxax, nfam);
               nfam = maxax;
            end
            resize = ternaryOp( nfam>6, 2, 1 ); % lots in family -> resize plots
            nr = resize*nap; nc = ceil(nfam/resize);
            
            max_count = 0;
            for jj=1:nfam
               numspikes = length( tseries.APstimes{ii}{jj} ); 
               peakinds  = round( tseries.APstimes{ii}{jj} / dt );
               peaks     = tseries.data{ii}(peakinds, jj);
               allpeaks  = [ allpeaks; peaks(:) ];
               alltimes  = [ alltimes; peakinds(:)*dt ];
               if min(peaks) < mv, mv = min(peaks); end
               if max(peaks) > Mv, Mv = max(peaks); end
               
               str = sprintf('AP %d, axon %d: %d spikes', ii, jj, numspikes);
               figure(ft);
               subplot(nr, nc, (ii-1)*nfam + jj)
                  plot( peakinds*dt, peaks, '.' );
                  title(str);
                  xlabel('Time (s)');
                  ylabel('Spike amp');
               % bin size lower bound is 20, & upper bound is 50
               nbins = min( [ 50 max([ ceil(numspikes/5) 20])] );
               figure(fh);
               subplot(nap, nfam, (ii-1)*nfam + jj)
                  [count, bins] = hist( peaks, nbins );
                  max_count = max(max_count, max(count)/sum(count));
                  bar( bins, count/sum(count) );
                  title(str);
                  xlabel('Spike amp');
                  ylabel('Probability');
            end
            
            % Adjust y-axis
            for jj = 1:nfam
               subplot(nap, nfam, (ii-1)*nfam + jj);
               ylim([0 max_count]);
            end
            
            ah = findobj(ft, 'type', 'axes');
            set(ah, 'ylim', [0 Mv], 'xlim', [1 nT]*dt);
            
         end
         ax = findobj(ft, 'type', 'axes');
         set(ax, 'ylim', [mv Mv]);
         ax = findobj(fh, 'type', 'axes');
         set(ax, 'xlim', [mv Mv]);
         
         fopts = {'fontweight', 'bold', 'fontsize', 16};
         set(ft, 'Visible', 'on');
         set(fh, 'Visible', 'on');
         figure; 
         subplot(2,1,1),
            plot( alltimes, allpeaks, 'o');
            title('Amplitude change for all spikes', fopts{:});
            xlim([1 nT]*dt);
            ylim([0 Mv]);
         subplot(2,1,2),
            nbins = min( [ 100 max([ ceil(length(allpeaks)/5) 20])] );
            [count, bins] = hist(allpeaks, nbins);
            bar(bins, count/sum(count));
            title('Amplitude probability for all spikes', fopts{:});
   end

end