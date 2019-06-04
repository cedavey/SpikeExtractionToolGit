% str = getCatchMEstring( ME, msg, noprint )
%
% Extract information from catch ME, & append to input message. Outputs
% string to stnd out in red text, unless noprint is set to true. 
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