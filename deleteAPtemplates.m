% tseries = deleteAPtemplates( tseries, method, method_params )
%
% Created by Artemio - 19/June/2019
function tseries = deleteAPtemplates( tseries, method, method_params )
   deleteAPs  = tseries.deleteAPs;
   
   % remove AP families that were merged
   tseries.APfamily( deleteAPs ) = [];
   tseries.data( :, deleteAPs ) = [];
   tseries = rmfield( tseries, 'deleteAPs' );
end