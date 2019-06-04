% [spikes,stimes] = runSpikeOperations(tseries, method, method_params)
function [spikes,stimes] = runSpikeOperations(tseries, method, method_params)

   switch lower(method)
      case 'merge spikes'
         bytemp = method_params.within_template.value;
         nT     = size(tseries.data{1}, 1); % get number of samples of 1st template
         nap    = length(tseries.data);

         % if merging by template, get a timeseries for each AP template,
         % else get a single timeseries for the entire dataset
         if bytemp
            % create a cell for each template
            spikes    = cell(nap, 1);
            stimes    = cell(nap, 1);
         else
            % create 1 template, and 1 family within the template
            spikes    = cell(1, 1);
            spikes{1} = zeros(nT, 1); 
            stimes    = [];
         end
         for ii=1:nap
            % if merging to AP template sum spikes across all families 
            % within each template, else sum across all families and then 
            % across all templates
            if bytemp
               spikes{ii} = sum( tseries.data{ii}, 2 ); 
               stimes{ii}{1} = sort( cell2mat( tseries.APstimes{ii}(:)' ) )';
               
            % merging all families and all templates
            else
               spikes{1}  = spikes{1} + sum( tseries.data{ii}, 2 );
               nspikes    = sum( cellfun( @length, tseries.APstimes{ii} ) );
               stimes(end+1 : end+nspikes) = cell2mat( tseries.APstimes{ii}(:)' );
            end
         end

         % if merging across all templates, sort spike times
         if ~bytemp
            % spike functions expect data to have structure of family
            % within each template, so create 1 template, with 1 family
            tmp = { sort( stimes ) }; clear stimes;
            stimes{1}{1} = tmp;
         end
   end

end