% RUNTIMEERRORHANDLER Deals with catched errors and logs them into a file 
% (different to the diary 'log_all.log', which saves everything that's 
% displayed on the command window).
% 
% Syntax: runtimeErrorHandler(MatlabError, action<optional>,message<optional>);
%
% Inputs:
%     MatlabError - catched error (ME)
%     action      - String, can be:
%                    'rethrow'(default) Logs the error and throws it again.
%                    'ignore'  Logs the error and then does nothing.
%                    'message' Logs the error and prints a message on the
%                              command window.
%     message     - If action is 'message', it prints this message to cmd.
%
% See also: try, catch, rethrow, error
%
% Created by Artemio - 19/June/2019
function runtimeErrorHandler(varargin)
% Will deal with the catched errors and log them to a file (different to
% the diary 'log_all.log', which saves everything that's displayed on the
% command window).
   if nargin >0, ME = varargin{1}; else, error('Not enough input arguments.');end
   if nargin > 1, action = varargin{2}; else, action = 'rethrow';end
   
   % Get location of log files
   path = getFilePath('log');

   %open file
   fid = fopen([path 'logFile.log'],'a+');
   % write the error to file
   % first line: message
   fprintf(fid,'\nError handled (%s) %s:\n%s\n',action,string(datetime),ME.message);
   % second line: matlab version
   mlversion = version;
   fprintf(fid,'MATLAB v%s\n',mlversion);
   % following lines: stack
   for e=1:length(ME.stack)
      fprintf(fid,' in %s at %i\n',ME.stack(e).name,ME.stack(e).line);
   end
   % close file
   fclose(fid);
   
   % Deal with the error
   switch action
      case 'rethrow'
         rethrow(ME);
         
      case 'ignore'
         return
         
      case 'message'
         if nargin < 3
            str = sprintf('\tFollowing error catched:\n\t%s\n\tFunction: %s > Line: %d.\n',ME.message, ME.stack(1).name, ME.stack(1).line);
         else
            str = sprintf('%s',varargin{3});
         end
         printMessage('on','Errors',str);
         
      otherwise
         error('I don''t know what to do with this error: %s', ME.message);
   end
end