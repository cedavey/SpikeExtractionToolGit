% handles = updateSETFigure(handles, tseries, tlim, vlim)
% Update figure panel on spike extraction tool (SET) to display current
% timeseries data
function handles = updateSETFigure(handles, tseries)
   % each voltage time series struct has the following fields available:
   %   tseries.name = name;
   %   tseries.type   = 'voltage';
   %   tseries.data   = voltage;
   %   tseries.time   = time;
   %   tseries.dt     = time(2) - time(1);
   %   tseries.params = [];

   if any(strcmpi(tseries.type, {'spike','rate'}))
      % plot spike & rate separately because data is in cells, rather than arrays
      handles = plotDataFromCell(handles, tseries);
   else
      handles = plotDataFromMatrix(handles, tseries);
   end



end

% Cell data has a timeseries for each cell
% Each AP template of a certain shape has a set of families with matching
% shape but different amplitudes. Plot these on the same axes.
function handles = plotDataFromCell(handles, tseries)
   tlim  = double(handles.data.tlim);
   vlim  = double(handles.data.vlim);
   panel = handles.plot_panel;
   %subplot(111); % resets axes
   arrayfun( @cla, findobj(handles.figure1,'type','axes' ) );

   max_ax= 3;      % max number of axes per figure - scrollbar to access others
   max_lines = 10; % max number of families per template (i.e. per axis)
   warn  = false;  % have given warning (remember so don't do lots of times for same thing)

   type  = tseries.type;
   data  = tseries.data;  % format: AP templates x 1
   Np    = length(data);  % num axes to include = number of AP templates
   if Np > max_ax % if too many axes to view add a scroll bar
      set(handles.scroll_axes_slider, 'Visible', 'on'); % no scrolling required
      scroll    = get(handles.scroll_axes_slider, 'Value');
      poss_ax   = (1:Np)';
      select_ax = poss_ax - round((Np - max_ax)*(1-scroll)); % scroll is max when want min axes so 1-scross
      plot_ax   = find( 1 <= select_ax & select_ax <= max_ax );
      Np        = max_ax;
      set(handles.scroll_axes_slider,'SliderStep',[min(1,2/length(data)) min(max_ax/length(data))]);
   else
      plot_ax = 1:Np;
      set(handles.scroll_axes_slider, 'Visible', 'off'); % no scrolling required
      set(handles.scroll_axes_slider, 'Value', 1); % reset scroll bar
   end
   nc  = 1;
   nr  = ceil(Np);
   dt  = tseries.dt;
   if tlim(2)<(tlim(1)+dt)
      tlim(2) = tlim(1) + dt;
   end
   if vlim(2) < (vlim(1)*1.1)
      vlim(2) = vlim(1)*2;
   elseif all(vlim==0)
      vlim = [0 1];
   end

   % I keep having float related problems so calc indices first
   sind = max( round(double(tlim(1)/dt)), 1 );
   eind = round( double(tlim(2)/dt) );
   time = toVec((sind:eind)*dt);

   % time & data keep being out by 1 sample, so just cheat...bad you!
   % - only do for arrays, not for event based spikes coz we'll extract
   % time for them later
   if isnumeric(data{1}) && length(time) > size(data{1},1)
      eind = size(data{1},1);
      time = time(1:eind);
   % can't do the else part coz we're going to set data limits to match axes
   % limits below, in the loop, so we don't apply it to more datasets than
   % we're actually going to display
%    elseif length(time)<size(data{1},1)
%       time = [time; time(end)+dt];
%       eind = eind + 1;
   end
   [tscale, tlabel] = getUnitScale(tlim(2)*10, 's');
tscale=1; tlabel='s';

   fopts= {'fontsize', 10, 'fontweight', 'bold'};
   for i=1:Np
      if ~isempty( data{ plot_ax(i) } )
         ax1 = subplot(nr, nc, i, 'Parent', panel); hold off;

         % spike data can be either in a timeseries, or in event based spikes
         if ~iscell( data{ plot_ax(i) } )
            x = data{ plot_ax(i) }(sind:eind, :);
         else
            % extract spikes for each axon family
            nfam = length( tseries.data{i} );
            x    = cell( 1, nfam );
            t    = cell( 1, nfam );
            for ff=1:nfam
               fam    = tseries.data{i}{ff};
               stimes = tseries.APstimes{i}{ff};
               valid  = time(sind) <= stimes & stimes <= time(eind);
               nV     = sum( valid );
               % to avoid a line joining one spike to the next making
               % random diagonals on the figure, add a zero at the beginning
               % & end of each spike, so the line is on the x-axis
               % - this needs a repeat of 1st & last time for each spike
               sptime = fam.stimes(:, valid);
               sp     = fam.spikes(:, valid);
               sptime = [ sptime(1,:); sptime; sptime(end,:) ];
               sp     = [ zeros(1,nV); sp; zeros(1,nV) ];
               t{ff}  = toVec( sptime );
               x{ff}  = toVec( sp );
            end
         end
         % only show up to max_lines different spikes with the mean AP
         % template (there can be 300 K spikes for a single template, which
         % destroys matlab). Randomise the lines else they can look really
         % wrong for the template --> the mean is not at all stationary
         numaxons = size( x, 2 ); % number of axons contributing to template
         if numaxons <= max_lines
            lind = 1:numaxons;
         else
            % get the first max_lines indices from numspikes
            lind  = randperm( numaxons, max_lines );
            if isnumeric(x)
               [~, idx] = sort(max(abs(x(:,lind))),'descend');
            else
               xMax = zeros(1,(numel(lind)));
               for ii = 1:numel(lind)
                  xMax(ii) = max(abs(cell2mat(x(lind(ii)))));
               end
               [~, idx] = sort(xMax,'descend');
            end
            lind = lind(idx);
         end

         [vscale, vlabel] = getUnitScale( vlim(2)*10, 'V' );

         numlines = min([max_lines, size(x,2)]);

         if isnumeric(x)
            try
               max_pks = max(x(:,lind));
               [~, sorted] = sort(max_pks, 'descend');
               lh = plot(ax1, time/tscale, x(:,sorted)/vscale);
            catch ME
               % keep having trouble with floats making time vector the wrong size
               lh = plot(ax1, (1:size(x,1))*dt/tscale, x(:,lind)/vscale);
               if ~warn % don't print over and over
                  str = sprintf('Warning: size of time and data didn''t match, reconstructing time vector\n');
                  printMessage('off','Comments*', str);
                  warn= true;
               end
            end
         % event-based plotting -> plot each axon family separately
         else
            cols = getColourMatrix( numlines );
            for ff=1:numlines
               plot(ax1, t{lind(ff)}, x{lind(ff)}, '-', 'color', cols(ff,:));
               hold on;
            end
         end
         xlim( tlim/tscale );
         ylim( vlim/vscale );
         % Count total found spikes to display the value on the title
         method_idx = ternaryOp(length(data) == size(tseries.APstimes,2), 2, 1); % 2 for moving average, 1 for k-means
         num_templates = length(data);% size(tseries.APstimes,method_idx);
         total_found_spikes = zeros(num_templates,1);
         try
             for j = 1:num_templates % debugging, not sure if num_templates is correct under every possibility
                 if method_idx > 1
                     for k = 1:numel(tseries.APstimes{1,j})
                        total_found_spikes(j) = total_found_spikes(j) + size(tseries.APstimes{1,j}{k},method_idx);
                     end
                 else
                     for k = 1:numel(tseries.APstimes{j})
                        total_found_spikes(j) = total_found_spikes(j) + size(tseries.APstimes{j}{k},method_idx);
                     end
                 end
             end
             set( get(ax1,'title'), 'String', sprintf('AP %d (%d fams, N = %d, Total = %d)', ...
                                                   plot_ax(i), size(x,2), total_found_spikes(plot_ax(i)), ...
                                                   sum(total_found_spikes)), fopts{:});
         catch E
               str = 'Couldn''t display the spike count here. Check in the workspace.\n';
               runtimeErrorHandler(E,'message',str);           
         end
         
         if i==1
            set( get(ax1,'Xlabel'), 'String', sprintf('Time (%s)', tlabel), fopts{:});
            switch tseries.type
               case 'spike'
                  set( get(ax1,'Ylabel'), 'String', sprintf('Voltage (%s)', vlabel), fopts{:} );
               case 'rate'
                  set( get(ax1,'Ylabel'), 'String', sprintf('Rate (%s)', vlabel), fopts{:} );
               otherwise
                  str = sprintf('Unknown time series type (%s) in updateSETfigure', tseries.type);
                  displayErrorMsg(str);
                  return;
            end
         end % end axis labelling if i==1
      end % end if data is not empty
   end % end for each plot_ax indices
end % end plotDataFromCell


function handles = plotDataFromMatrix(handles, tseries)
   tlim = double(handles.data.tlim);
   vlim = double(handles.data.vlim);
   panel = handles.plot_panel;
   % subplot(1,1,1); % resets axes
   arrayfun(@cla,findobj(handles.figure1,'type','axes'))

   max_ax= 9;
   warn  = false;   % have given warning (remember so don't do lots of times for same thing)
   max_lines = 500; % max number of spikes contributing to mean AP template

   type  = tseries.type;
   data  = tseries.data;   % format: time x dimension
   
   % extended code to allow templates to have diff lengths 
   if iscell( data )
      uniqueAPLength = true;
      Np = length(data);  % num axes to include
      nT = cellfun( @(d) size( d, 1 ), data );
   else
      uniqueAPLength = false;
      [nT, Np] = size( data );  % num axes to include
   end
      
   if Np > max_ax % if too many axes to view add a scroll bar
      set(handles.scroll_axes_slider, 'Visible', 'on'); % no scrolling required
      scroll    = get(handles.scroll_axes_slider, 'Value');
      poss_ax   = (1:Np)';
      select_ax = poss_ax - round((Np - max_ax)*(1-scroll)); % scroll is max when want min axes so 1-scross
      plot_ax   = find( 1 <= select_ax & select_ax <= max_ax );
      Np        = max_ax;
      fontsize  = 10;
   else
      plot_ax = 1:Np;
      set(handles.scroll_axes_slider, 'Visible', 'off'); % no scrolling required
      set(handles.scroll_axes_slider, 'Value', 1); % reset scroll bar
      fontsize  = 12;
   end
   nc = ceil(sqrt(Np));
   nr = ceil(Np/nc);
   if strcmp('gui', handles.f.uiType)
      set(handles.scroll_axes_slider,'SliderStep',[min(1,(1+max_ax)/Np) min(1,(1+Np*max_ax)/Np)]);
   end
   dt = tseries.dt;
   if tlim(2)<(tlim(1)+dt)
      tlim(2) = tlim(1) + dt;
   end
   if vlim(2)<(vlim(1)+dt)
      vlim(2) = vlim(1)*2;
   end

   % i keep having float related problems so calc indices first
   sind = max( round( double( tlim(1) / dt ) ), 1);
   eind = round ( double( tlim(2) / dt) );
   time = toVec( (sind:eind) * dt );
   % time & data keep being out by 1 samples, so just cheat...bad you!
   if eind > max( nT )
      time = time( 1:end - ( eind-size(data,1) ) );
      eind = eind - ( eind - size(data,1 ));
   end
   if uniqueAPLength
%       data = cellfun( @(d) d(sind:eind, :), data );
      [ vscale, vlabel ] = getUnitScale( max( cellfun( @max, data ) )*10, 'V' );
   else
      data = data(sind:eind, :);
      [ vscale, vlabel ] = getUnitScale( max(data(:))*10, 'V' );
   end
   if length(time) < size(data,1)
      time = [time; time(end)+dt];
      eind = eind + 1;
   end
   [ tscale, tlabel ] = getUnitScale( tlim(2)*10, 's' );
tscale=1; tlabel='s';

   fopts= {'fontsize', fontsize, 'fontweight', 'bold'};
   figure(handles.figure1); % Change focus to main window
   if strcmp('app', handles.f.uiType)
      panel.AutoResizeChildren = 'off';
      panel.Children(4).delete;
   else
      subplot(1,1,1, 'Parent', panel); % Clear the axes to prevent having empty little subplots
   end
   for i=1:Np
      ax1 = subplot(nr, nc, i, 'Parent', panel); hold off;
      if strcmp('app', handles.f.uiType)
         try
            close(1);
         catch
            nan;
         end
      end % In app, there's a new Figure window opening for some reason.
      % try plotting data at chosen scale
      try
         if uniqueAPLength
            lh = plot( ax1, time(1:nT(plot_ax(i)))/tscale, data{plot_ax(i)}/vscale );
            numT = nT(plot_ax(i));
         else
            lh = plot( ax1, time/tscale, data(:, plot_ax(i))/vscale );
            numT = nT;
         end
      catch ME
         % keep having trouble with floats making time vector the wrong size
         if uniqueAPLength
            lh = plot( ax1, (1:size(data{plot_ax(i)}))*dt/tscale, data{plot_ax(i)}/vscale );
         else
            lh = plot( ax1, (1:size(data(:, plot_ax(i))))*dt/tscale, data(:, plot_ax(i))/vscale );
         end
         if ~warn % don't print over and over
            str = sprintf('Warning: size of time and data didn''t match, reconstructing time vector\n');
            printMessage('off','Comments*', str);
            warn= true;
         end
      end
      set( ax1, 'xlim', tlim/tscale );
      set( ax1, 'ylim', vlim/vscale );

      set( get(ax1,'Xlabel'), 'String', sprintf('Time (%s)',    tlabel), fopts{:} );
      set( get(ax1,'Ylabel'), 'String', sprintf('Voltage (%s)', vlabel), fopts{:} );

       switch tseries.type
          case 'voltage'
             set( get(ax1,'Xlabel'), 'String', sprintf('Time (%s)',    tlabel), fopts{:});
             set( get(ax1,'Ylabel'), 'String', sprintf('Voltage (%s)', vlabel), fopts{:});

          case 'ap'
             fontsize  = 10; % make font smaller cuz looks silly!
             set(lh, 'color', 'k', 'linewidth', 4);
             try
                hold on;
                % nlines = min( max_lines, size( tseries.APfamily{plot_ax(i)}, 2) );
                size_lines = size( tseries.APfamily{ plot_ax(i)}, 2 );
                nlines     = min( max_lines, size_lines );
                % Show the last max_lines lines, instead of the first ones
                lind       = max( 1, (size_lines - max_lines) ) : size_lines;

                plot( ax1, time(sind:min(eind,numT))/tscale, tseries.APfamily{plot_ax(i)}(sind:min(eind,numT),lind)/vscale, 'linewidth', 1 );
                vlim(2)    = max( toVec( tseries.APfamily{plot_ax(i)}(sind:min(eind,numT),lind) ) );
                % if APs normalised we know they'll be btwn -3/3 so make lims the same
                if tseries.params.normalise_aps.value
                   ylim( [-3 3] );
                else
                   ylim( vlim/vscale );
                end
                hold off;
                uistack(lh, 'top'); % put mean on top
    %                chH = get(ax1,'Children');
    %                set(gca,'Children',[chH(end);chH(1:end-1)]);
                N = size( tseries.APfamily{plot_ax(i)}, 2 );
                set( get(ax1,'Title'), 'String',  sprintf('AP %d (N=%d)', plot_ax(i), N), 'fontweight','bold','fontsize',10 );
                set( get(ax1,'Ylabel'), 'String', sprintf('Voltage (%s)', vlabel), 'fontweight','normal','fontsize',10 );

                if tseries.params.normalise_aps.value
                   set( get(ax1,'Ylabel'), 'String', sprintf('Normalised'), 'fontweight','normal','fontsize',10 );
                end

             catch ME
                if strcmp('Invalid or deleted object.', ME.message)
                   set(handles.figure1,'CurrentAxes',ax1)
                   if ~ishold
                      hold(ax1,'on');
                   end
                   plot((1:size(data(:,...
                      plot_ax(i))))*dt/tscale, data(:, plot_ax(i))/vscale,...
                      'k','LineWidth',3);
                else
                   runtimeErrorHandler(ME,'message');
                end
            end %try
        end %switch
   end %for
end %function
