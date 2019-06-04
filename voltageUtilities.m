function new_tseries = voltageUtilities(tseries, method, params)
   dt  = tseries.dt;

   switch lower( method )
      case 'downsample'
         factor = params.factor.value; 
         switch lower( params.method.value )
            case 'ditch'
               voltage = tseries.data(1:factor:end); 
            case 'average'
               meanfn = @(block) mean2(block.data(:));
               voltage = blockproc( tseries.data(:), [factor 1], meanfn );
         end

         % make new tseries structure with voltage, time, dt
         new_tseries.data   = voltage; 
         new_tseries.dt     = dt * factor; 
         new_tseries.time   = tseries.time(1) : new_tseries.dt : tseries.time(end);
         new_tseries.params = params;
       
      case 'truncate'
         type   = params.type.value; % samples or seconds
         amount = params.amount.value; % num sample or seconds to ditch
         if strcmpi( type, 'seconds' )
            amount = round( amount / dt ); % convert seconds to num samples
         end
         if amount > length( tseries.data )
            displayErrorMsg( 'Cannot truncate by more than length of data, ignoring request' );
            new_tseries = []; 
            return;
         end
         switch lower( params.method.value )
            case 'start'
               voltage = tseries.data(amount+1:end);
               new_tseries.time = tseries.time(amount+1:end);
            case 'end'
               voltage = tseries.data(1:end-amount);
               new_tseries.time = tseries.time(1:end-amount);
         end

         % make new tseries structure with voltage, time, dt
         new_tseries.data   = voltage; 
         new_tseries.dt     = dt; 
         new_tseries.params = params;

   end % end utilities method

end















