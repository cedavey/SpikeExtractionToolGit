% bool = iseven(x)
% Returns true if the length of x is even. You can input the length or the
% actual vector itself.
function bool = iseven(x)
   if ~isscalar(x)
      x = length(x);
   end
   if x/2==round(x/2)
      bool = true;
   else
      bool = false;
   end
end