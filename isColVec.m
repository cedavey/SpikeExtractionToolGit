% y = isColVec(v)
function y = isColVec(v)
   if nargin==0
      help isColVec; 
      return;
   end
    if ~any( size(v) == 1 )
%         cprintf('Error', 'Input is not a vector at all');
        y = false;
        return;
    end

    if size(v,1)>size(v,2)
        y = true;
    else
        y = false;
    end
end
