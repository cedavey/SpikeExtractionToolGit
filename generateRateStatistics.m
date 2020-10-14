% generateRateStatistics(tseries, method, method_params)
% Generate statistics on rate timeseries for each AP family. Currently
% the only statistic is autocorrelation. 
function generateRateStatistics(tseries, method, method_params)
   maxlag  = method_params.max_lag.value;
   dt      = tseries.dt;
   numlags = ceil(maxlag/dt); % lags given in ms
   numlags = ternaryOp( numlags<=1, maxlag, numlags );
   maxcols = 10; 
   
   if numlags > 100 % if lots of lags check there wasn't a mistake
      str = sprintf('Generating autocorrelation of length %d, continue?', round(maxlag/dt));
      response = userConfirmation(str,'Number of lags');
      switch lower(response)
        case 'yes'
           % do nothing
        case 'no'
           return;
      end
   end

   % for each AP templates, plot the autocorrelation for each family 
   nr = length(tseries.data);
   fh = figure; set(fh, 'Visible', 'off');
   
   for ii=1:nr
      APrates  = tseries.data{ii}; 
      [nT, nc] = size(tseries.data{ii}); 
      nc = min( nc, maxcols );
      
      for jj=1:nc
         % only include when spiking has started, up until it's finished
         x    = APrates(:, jj); 
         sind = find( x > 0, 1, 'first' );
         eind = find( x > 0, 1, 'last' );
         x    = x( sind:eind );

         try
            [rho, lags] = autocorr(APrates(:, jj), numlags);
         catch
            [rho, lags] = xcorr(APrates(:, jj), numlags);
         end
         str = sprintf('AP %d (%d)', ii, jj);
         subplot(nr, nc, (ii-1)*nc + jj)
%             bar(lags(end-numlags:end)*dt, rho(end-numlags:end));
            stem( lags*dt, rho );
            title( str );
            ylabel( 'Correlation' );
            xlabel( 'Time (s)' );
            xlim([-1 numlags]);
      end
      
      set( gcf, 'name', sprintf( 'Rate autocorrelation (max %d families per template shown)', maxcols ) );
   end
   set(fh, 'Visible', 'on');
end