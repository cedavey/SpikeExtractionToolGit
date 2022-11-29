% generateSpikeStatistics(tseries, method, method_params)
% Generate statistics on rate timeseries for each AP family. Currently
% the only statistic is autocorrelation. 
function generateSpikeStatistics(tseries, method, method_params)
   nap     = length( tseries.data );
   dt      = tseries.dt;
   
   try
      % determine max number of rows & columns, and actual number of rows &
      % columns, when plotting each stat across templates & families
      maxcols = 6; % max plot columns per template
      maxrows = 3; % max plot rows per template
      nfam    = cellfun( @(d) size( d, 2 ), tseries.data );
      maxrows = min( ceil( max(nfam) / maxcols ), maxrows ); % might not need all rows 
      maxax   = maxcols * maxrows; % max axons within each family to display
      if any( nfam > maxrows * maxcols )
         str  = sprintf('\tOnly displaying a max of %d families from each template \n', maxax );
         cprintf( 'Keywords*', str );
      end
      nfam( nfam > maxax ) = maxax;
      % calc num cols & num rows in each figure - same for each statistic so calculate once
      ncols   = min( max( nfam ), maxcols );    % e.g. max of 2 rows per AP family, max 6 cols per row
      AProws  = ceil( nfam / maxcols );         % number of rows allocated to this template
      nrows   = sum( AProws );                  % sum of rows for each template

      % the colour matrix for each template is the same here as it is for the
      % debugging plots of spike peak mean & variance for
      % extractSpikeUsingTemplate, to allow for visual comparison
      switch lower(method)
         %% spike raster
         case 'raster'
            fopts = {'fontweight', 'bold', 'fontsize', 16};
            fh    = figure; set( fh, 'Visible', 'off' ); % figure for histogram of amplitudes

            % for each AP templates, plot the autocorrelation for each family 
            % cols  = getColourMatrix(nap); % get each template a unique colour
            index = 0;
            for ii=1:nap
               cols = getColourMatrix( nfam(ii) );

               for jj=1:nfam(ii)
                  stimes = round( tseries.APstimes{ii}{jj} );
                  index  = index + 1;
                  figure( fh );
                     hold on;
                     plot(stimes, index * ones(size(stimes)), '.', 'color', cols(jj,:) );
               end
            end
            set(fh, 'Visible', 'on'); % figure for histogram of amplitudes
   %          set(gca, 'xcolor', [1 1 1]);
            set(gca, 'ycolor', [1 1 1]);
            xlabel('Time (seconds)', fopts{:});
            xlim( [tseries.time(1) tseries.time(end)] )
            ylim( [1 index] );
            title( sprintf( 'Raster: AP %d axon %d', ii, jj ) ) ;

         case 'psth'
            fopts = { 'fontweight', 'bold', 'fontsize', 12 };
            fh = figure; set( fh, 'Visible', 'off' ); % figure for histogram of amplitudes
            numsecs  = method_params.time_post_spike.value;
            num_bins = method_params.num_bins.value;
            numsamp  = round( numsecs / dt );
            dbin     = numsecs / num_bins;
            timebins = dbin:dbin:(numsecs+dbin); % add list bin to catch everything above

            % for each AP templates, plot the autocorrelation for each family 
            for ii=1:nap
               cols = getColourMatrix( nfam(ii) );

               % start new APs axes at beginning of new row
               fi = ternaryOp( ii==1, 1, sum( AProws(1:ii-1) ) * ncols + 1 ); 
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
                     subplot(nrows, ncols, fi), fi=fi+1;
      %                   plot(timebins(1:end-1), timecount, '.', 'color', cols(ii,:) );
                        b = bar( timebins(1:end-1), timecount, 'facecolor', 'flat' ); % weird to control colour
                        b.FaceColor = cols(jj,:);
                        set(get(gca,'children'),'facecolor',cols(jj,:));
                        box off;
                        xlim( [ timebins(1)-(timebins(2)-timebins(1)) timebins(end) ] );
                        xlabel('Time (s)', fopts{:});
                        xlabel('P(spike)', fopts{:});
                        title( sprintf('PSTH: AP %d axon %d', ii,jj) ) ;
               end
            end
            set(fh, 'Visible', 'on'); % figure for histogram of amplitudes

         %% ISI
         case 'interspike interval'
            % for each AP templates, plot the ISI for each family          
            ft = figure; set( ft, 'Visible', 'off' );
            num_bins = method_params.num_bins.value;
            max_isi  = method_params.max_isi.value;

            for ii=1:nap
               % start new APs axes at beginning of new row          
               fi      = ternaryOp( ii==1, 1, sum( AProws(1:ii-1) ) * ncols + 1 ); % ncols instead of maxcols
               max_cnt = 0; ah = zeros( nfam(ii), 1 );
               cols    = getColourMatrix( nfam(ii) );
               for jj=1:nfam(ii)
                  stimes = tseries.APstimes{ii}{jj};
                  isi    = diff(stimes);
                  isi(isi>max_isi) = [];
                  [count, bins] = hist(isi, num_bins);
                  % figure(fa); subplot(nr, nc, fi); autocorr(isi,'NumLags',size(isi,2)-1);
                  count   = count / sum(count);
                  max_cnt = ternaryOp( max(count) > max_cnt, max(count), max_cnt );
                  str     = sprintf('AP %d, axon %d: %d spikes', ii, jj, length(stimes));
                  figure(ft);
                     ah(jj) = subplot(nrows, ncols, fi); fi=fi+1;
                        b = bar( bins, count, 'facecolor', 'flat' ); % weird to enable controlling color
                        b.FaceColor = cols(jj,:);
                        % set(get(gca,'children'),'facecolor',cols(ii,:));
                        title(str);
                        xlabel('Time (s)');
                        ylabel('Probability');
                        title( sprintf('ISI: AP %d axon %d', ii,jj) ) ;
                        xlim([0 max_isi]); % To show the bins from 0 to the user chosen upper limit.
               end
               set( ah, 'ylim', [0 max_cnt] ); 
            end
            set(ft, 'Visible', 'on');

         %% amplitude change over time
         case 'amplitude change'
            % for each AP templates, plot the autocorrelation for each family 
            ft    = figure; set(ft, 'Visible', 'off'); % figure for amp change over time
            fh    = figure; set(fh, 'Visible', 'off'); % figure for histogram of amplitudes
            mv    =  Inf;
            Mv    = -Inf;
            vmin  = min( cellfun( @(v) min( v(:) ), tseries.data ) ) * 0.8;
            vmax  = max( cellfun( @(v) max( v(:) ), tseries.data ) ) * 0.8;
            nbins = 40;
            dbin  = ( vmax - vmin ) / nbins;
            bins  = vmin : dbin : vmax; 

            allpeaks = []; alltimes = [];
            for ii=1:nap
               nT      = size( tseries.data{ii}, 1 ); 
               max_cnt = 0; ah1 = zeros( nfam(ii), 1 ); ah2 = ah1;
               fi      = ternaryOp( ii==1, 1, sum( AProws(1:ii-1) ) * ncols + 1 );            
               for jj=1:nfam(ii)
                  cols      = getColourMatrix( nfam(ii) );
                  numspikes = length( tseries.APstimes{ii}{jj} ); 
                  peakinds  = round(  tseries.APstimes{ii}{jj} / dt );
                  peaks     = tseries.data{ii}(peakinds, jj);
                  allpeaks  = [ allpeaks; peaks(:) ];
                  alltimes  = [ alltimes; peakinds(:)*dt ];
                  if min(peaks) < mv, mv = min(peaks); end
                  if max(peaks) > Mv, Mv = max(peaks); end

                  str = sprintf('AP %d, axon %d: %d spikes', ii, jj, numspikes);
                  figure(ft);
                     ah1(jj) = subplot(nrows, ncols, fi ); 
                        plot( peakinds*dt, peaks, '.', 'Color', cols(jj,:) );
                        title(str);
                        xlabel('Time (s)');
                        ylabel('Spike amp');
                        xlim( [ 0 nT*dt ] );

                  figure(fh);
                     ah2(jj) = subplot( nrows, ncols, fi ); fi=fi+1;
                        count       = hist( peaks, bins );
                        max_cnt     = max( max_cnt, max(count)/sum(count) );
                        b           = bar( bins, count/sum(count), 'facecolor', 'flat' ); % weird to enable controlling color
                        b.FaceColor = cols(jj,:);
                        title(str);
                        xlabel('Spike amp');
                        ylabel('Probability');
                        title( sprintf('Amp: AP %d axon %d', ii,jj) ) ;
                        axis square;
               end
               set( ah2, 'ylim', [0 max_cnt] ); 

               % ah = findobj( ft, 'type', 'axes' );
               % set( ah, 'ylim', [0 Mv], 'xlim', [1 nT]*dt );

            end % end for each AP template
            ax = findobj( ft, 'type', 'axes' );
            set(ax, 'ylim', [mv Mv]*0.8);
            
            % ax = findobj( fh, 'type', 'axes' );
            % set(ax, 'xlim', [mv Mv]);

            fopts = { 'fontweight', 'bold', 'fontsize', 16 };
            set( ft, 'Visible', 'on' );
            set( fh, 'Visible', 'on' );
            figure; 
               subplot(211),
                  plot( alltimes, allpeaks, '.');
                  title( 'Amplitude of all spikes', fopts{:} );
                  xlim( [1 nT]*dt );
                  ylim( [0 Mv] );
               subplot(212),
                  count = hist( allpeaks, bins );
                  bar( bins, count/sum(count) );
                  title( 'Amplitude pdf for all spikes', fopts{:} );
      end

   catch ME
      str = getCatchMEstring( ME );
      displayErrorMsg( str );
   end
end