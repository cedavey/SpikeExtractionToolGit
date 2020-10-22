% tseries = deleteAPtemplates( tseries, method, method_params )
%
% Created by Artemio - 19/June/2019
function tseries = deleteAPtemplates( tseries, method, method_params )
   try
      deleteAPs  = tseries.deleteAPs;   
      % remove AP families that were merged
      tseries.APfamily( deleteAPs )   = [];
      if iscell( tseries.data )
         tseries.data( deleteAPs )    = [];
         tseries.time( deleteAPs )    = [];
      else
         tseries.data( :, deleteAPs ) = [];
      end
      tseries = rmfield( tseries, 'deleteAPs' );
      
   catch ME
      if strcmp('Reference to non-existent field ''deleteAPs''.',ME.message)
         displayErrorMsg('The template you are trying to delete does not exist.');
         tseries = [];
      else
         rethrow(ME);
      end
   end
end