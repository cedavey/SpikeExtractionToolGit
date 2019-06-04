% x = toVec(x)
% Forces x to column vector, regardless of shape
function x = toVec(x)
   if iscell(x)
      x = x{:};
   else
      x = x(:);
   end
end
