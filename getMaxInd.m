% ind = getMaxInd( x, dim, offset )
% 
% Return the index of the maximum value in x - allows inline equations
% 
% Inputs:
%      x - data to find max for, can be numeric or cell
%    dim - dimension to find max along (defaults to 1)
% offset - can offset along dimension when finding max, but still returns
%          index relative to 1 (defaults to no offset)
%          (supports up to 3 dimensions)
function ind = getMaxInd( x, dim, offset )
   if nargin<2, dim = 1; end
   if nargin<3, offset = 0; end
   
   if iscell(x) 
      dimsz = isDim( x{1} );
   else
      dimsz = isDim( x );
   end
   if isnumeric( x )
      switch dimsz
         case 1
            ind = getInd( x(1+offset:end-offset,:,:), dim );
         case 2
            ind = getInd( x(:,1+offset:end-offset,:), dim );
         case 3
            ind = getInd( x(:,:,1+offset:end-offset,:), dim );
      end
   elseif iscell( x )
      switch dimsz
         case 1
            ind = cellfun( @( c ) getInd( c(1+offset:end-offset,:,:), dim ), x );
         case 2
            ind = cellfun( @( c ) getInd( c(:,1+offset:end-offset,:), dim ), x );
         case 3
            ind = cellfun( @( c ) getInd( c(:,:,1+offset:end-offset), dim ), x );
      end
      
   end
   ind = ind + offset;
end

function ind = getInd( x, dim )
   [~, ind] = max( x, [], dim );
end