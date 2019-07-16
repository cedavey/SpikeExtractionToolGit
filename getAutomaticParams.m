% GETAUTOMATICPARAMS Returns the ideal parameters for the chosen method
% that fits the input voltage.
%
% Syntax: params = getAutomaticParams(tseries, method);
%
% Inputs:
%     tseries	- The voltage data
%     methods	- Method with which the data is being analyzed
%
% Output:
%     params	- Automatic parameters
%
% Artemio - 16 July 2019
%
function params = getAutomaticParams(tseries, method)
   % Check data type
   switch lower(tseries.type)
      case 'voltage'
         % Check method
         switch lower(method)
            case 'particle filter'
               
            otherwise
               % There is no automatic parameter configuration for the current
               % data type and method combination.
               str = sprintf('\tThere is no automatic parameter configuration for the data type: ''%s'' and method: ''%s''\n\tUsing the default parameters\n',tseries.type, method);
               printMessage('off','Errors',str);
         end
      otherwise
         % There is no automatic parameter configuration for the current
         % data type and method combination.
         str = sprintf('\tThere is no automatic parameter configuration for the data type: ''%s'' and method: ''%s''\n\tUsing the default parameters\n',tseries.type, method);
         printMessage('off','Errors',str);
   end
end