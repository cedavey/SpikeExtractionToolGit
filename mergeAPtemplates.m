% tseries = mergeAPtemplates( tseries, method, method_params )
function tseries = mergeAPtemplates( tseries, method, method_params )
   toid = method_params.template_to_merge_with.value;
   mergeAPs  = tseries.mergeAPs;
   
   for ii=1:length(mergeAPs)
      fromid = mergeAPs( ii );
      tseries.APfamily{ toid } = [tseries.APfamily{ toid } tseries.APfamily{ fromid }];
   end
   
   % recalculate mean AP for template that was merged in to
   tseries.data( :, toid ) = mean( tseries.APfamily{toid}, 2 );

   % remove AP families that were merged
   tseries.APfamily( mergeAPs ) = [];
   tseries.data( :, mergeAPs ) = [];
   tseries = rmfield( tseries, 'mergeAPs' );
end