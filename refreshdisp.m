% refreshdisp( str, prevstr, firstprint )
%
% Update the most recently output string. If it's the first time printing
% the string it's just printed as usual. 
function refreshdisp(str, prevstr, firstprint)

   if ~exist('firstprint','var') 
       firstprint = isempty(prevstr);
   end

   if firstprint
       fprintf(str)
   else
       fprintf( char( 8 * ones(1,length(prevstr))) );
       fprintf(str);
   end
   
end