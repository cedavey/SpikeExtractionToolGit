function printMessage( logMessage, varargin )
% PRINTMESSAGE displays styled formatted text in the Command Window
% 
% Syntax:
%    printMessage( logMessage, style, format, ... )
%
% Description:
% Displays a message in the command window using cprintf with the option of
% logging it to the log file.
%
% Inputs:
%   logMessage - if 'off', it doesn't save the printed message in the 
%                           log file. Anything else will save this
%                           message in the, log file.
%   all other inputs are passed to cprintf;  cprintf( format, str )
% 
% See also: cprintf
%
% Artemio - 12/June/2019

% If logMessage is false, turn off the log file.
isLogCurrentlyActive = strcmp('on',get(0,'Diary'));
if isLogCurrentlyActive && strcmp('off',logMessage), diary off; end

try
   % The string is formatted already, escaped characters, like %% now
   % appear like a single %. 
   format = varargin{1};
   str    = varargin{2};
   idx    = strfind(str,'%');
   if ~isempty(idx) && ~strcmp(str(idx + 1),'%')
      str(idx + 1 : end + 1) = ['%' str(idx+1:end)];
   end
   % Send the text and options to the command window
   cprintf(format, str);
catch ME
   % If an error, reactivate the diary
   if isLogCurrentlyActive, diary on; end
   rethrow(ME);
end

if isLogCurrentlyActive, diary on; end

end