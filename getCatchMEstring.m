% GETCATCHMESTRING - Extract information from catch ME, & append to input message.
% Outputs string to stnd out in red text, unless noprint is set to true. 
%
% str = getCatchMEstring( ME, msg, noprint )
%
% 
%
% Inputs:
%   ME      - catch exception
%   msg     - message to include in error output
%   noprint - if true don't print message to stdout
%
% Outputs:
%   str     - exception message str with line & function causing the error
function str = getCatchMEstring( ME, msg, noprint )
   if nargin<3 || isempty(noprint)
      noprint = false;
   end
   if nargin<2, msg = []; end
   if nargin<1 
      msg = 'Error';
   end
   str = sprintf( '%s: %s (line %d, function %s)\n', ...
                   msg, ME.message, ME.stack(1).line, ME.stack(1).name );
   if ~noprint
      cprintf( 'Errors*', str );
   end
end