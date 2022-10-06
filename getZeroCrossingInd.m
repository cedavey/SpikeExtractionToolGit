% ind = getZeroCrossingInd( x, dim, offset )
% 
% Return the index of the first zero crossing, from either positive to 
% negative or negative to positive - allows inline equations &
% cell inputs. 
% If multiple instances of max value the index of the first instance is
% returned. 
%
% x can be a matrix, in which case the zero crossing index will be returned
% for each vector, with dimension of interest specified by dim. If there is
% no zero crossing, ind = 0. 
% 
% Inputs:
%      x - data to find zero crossing for, can be numeric or cell
%    dim - dimension to find max along (defaults to 1)
% offset - can offset along dimension when finding zero crossing, but still 
%          returns index relative to 1 (defaults to no offset)
%          (supports up to 3 dimensions)
function ind = getZeroCrossingInd( x, dim, offset )
   if nargin<2,    dim = 1; end
   if nargin<3, offset = 0; end
   
   if iscell(x) 
      dimsz = isDim( x{1} );
   else
      dimsz = isDim( x );
   end
   if isvector( x )
      isrowvec = size( x, 1 ) < size( x, 2 );
      if isrowvec
         x = x';
      end
   end
   if isnumeric( x )
      switch dimsz
         case 1
            ind = getInd( x(1+offset:end-offset,:,:),   dim );
         case 2
            ind = getInd( x(:,1+offset:end-offset,:),   dim );
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
   if isvector( x )
      if isrowvec
         ind = ind';
      end
   end
end

function ind = getInd( x, dim )
   % if x is a vector we can simply use find directly
   if isvector(x)
      if x(1)>0
         ind = find( x < 0, 1, 'first' );
      else
         ind = find( x > 0, 1, 'first' );
      end
      if isempty(ind), ind = 0; end
      return;
   end

   % if x is an array we need to operate on vectors, with dimension of
   % interest specified by dim - use permute to get this dimension first
   Ndims = numel(size(x)); % total number of dimensions
   dims  = 1:Ndims;        % dimensions as a list
   dims(dim) = [];         % separate the dimension of interest from list
   X     = permute( x, [dim dims]);
   szX   = size(X); 
   Ncol  = szX(1);
   X     = reshape( X, Ncol, prod(szX(2:end)) );

   % doesn't return a result for vectors with no zero crossing
   ind   = table2array( rowfun( @(d) find(d<0,1,'first'), table(X) ) );

   ind = arrayfun( @(i) find(X(i,:)<0, 1, 'first'), 1:Ncol , 'UniformOutput', false, 'uni', false);
   ind = cellfun( @(i) ternaryOp( isempty(i), 0, i), ind );
end





