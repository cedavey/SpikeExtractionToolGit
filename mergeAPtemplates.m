% tseries = mergeAPtemplates( tseries, method, method_params )
function tseries = mergeAPtemplates( tseries, method, method_params )
   toid     = method_params.template_to_merge_with.value; % template to merge into
   mergeAPs = tseries.mergeAPs; % template(s) to merge from 
   
   % recalculate mean AP for template that was merged in to
   if iscell( tseries.data )
      % all template IDs
      allids = [toid; mergeAPs(:) ];
      % templates are allowed to have different lengths
      [ tseries.data{ toid }, tseries.APfamily{ toid } ] = ...
                      mergeTemplates( tseries.data( allids ), tseries.APfamily( allids ) );
      tseries.data( mergeAPs ) = [];
      tseries.time( mergeAPs ) = [];    
      % gotta update time since longest template may have been deleted
      if iscell( tseries.data )
         for ti=1:length(tseries.data)
            tseries.time{ti} = (1:size(tseries.data{ti},1))'*tseries.dt;
         end
      end
   else
      % templates are all forced to have the same length
      for ii=1:length(mergeAPs)
         fromid = mergeAPs( ii );
         tseries.APfamily{ toid } = [tseries.APfamily{ toid } tseries.APfamily{ fromid }];
      end
      tseries.data( :, toid ) = mean( tseries.APfamily{toid}, 2 );
      tseries.data( :, mergeAPs ) = [];
   end

   % remove AP families that were merged
   tseries.APfamily( mergeAPs ) = [];
   tseries = rmfield( tseries, 'mergeAPs' );
end

% merge AP templates that are highly correlated
% component APs for each template are already centred around the peak
function [mergedTemp, mergedComp] = mergeTemplates( templates, componentAPs )
   nTemp    = length( templates );
   peakind  = getMaxInd( templates ); % find index of template peak 

   nT       = cellfun(@length, templates );    % length of template
   tbefore  = max( peakind );
   tafter   = max( nT - peakind );
   
   if tbefore == 0 || tafter == 0
      str = sprintf( 'identifyUniqueTemplates:mergeTemplates:error: tbefore=%d, tafter=%d\n', tbefore, tafter );
   end
   
   % append NaNs to template spikes or new spike to make them the same size
   prependfn= @( sp, N ) [ NaN( N, size(sp,2)); sp ];
   appendfn = @( sp, N ) [ sp; NaN( N, size(sp,2)) ];
   
   try
      tempmean = zeros( tbefore + tafter, nTemp );
   catch ME
      str = sprintf( 'identifyUniqueTemplates:mergeTemplates:error: tbefore=%d, tafter=%d\n', tbefore, tafter );
      cprintf( 'Keywords*', str );
      runtimeErrorHandler( ME );
   end
   
   for i=1:nTemp
      % if this template has fewer samples before peak, prepend NaNs 
      if peakind(i) < tbefore  
         componentAPs{i} = prependfn( componentAPs{i}, tbefore - peakind(i) );
      end
      if nT(i) - peakind(i) < tafter
      % if this template has fewer samples after peak, append NaNs 
         componentAPs{i} = appendfn( componentAPs{i}, abs( nT(i) - peakind(i) - tafter ) );
      end
      tempmean(:,i) = nanmean( componentAPs{i}, 2 );
   end
   % mergedComp = horzcat( componentAPs{:} );
   mergedComp = [ componentAPs{:} ]; 
   
   % now take an ensemble average of template, where the number of samples
   % contributing to each time point may change across the spike, depending
   % on the length of the constituent spikes 
   mergedTemp = nanmean( mergedComp, 2 ); 
   
end

