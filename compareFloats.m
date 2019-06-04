% same = compareFloats(num1, num2, tol)
% Compare two numbers, including floats. Consider them equal if they're the
% same to within 'tol' tolerance. By default the difference between them is
% calculated using subtraction, but you can use the function handle input
% to change this to division, or whatever comparison technique is preferred
% Inputs:
%  num1 - number of any numerical datatype to compare (can be a vector)
%  num2 - number of any numerical datatype to compare (can be a vector)
%
%         If both inputs are vectors they must be the same size to allow
%         pairwise comparison. 
%
%  tol  - tolerance, within which to consider the numbers equal
%  fh   - function handle to allow different types of comparison. The
%         function must allow 2 inputs & return a single boolean.
%         Comparison is done in each direction, i.e. with each number
%         having a turn at being input to the function first.
%         Note: if fh=='percent' then function returns true if num1 &
%         num2 are within tol % of each other
% Outputs:
%  same - boolean, true if the numbers are the same
function same = compareFloats(num1, num2, tol, fh)
   if nargin==0, help compareFloats; return; end
   if nargin<3 || isempty(tol), tol=1e-3; end
   same = false;
   if (nargin<4 || isempty(fh)) && isscalar(num1) && isscalar(num2)
      same = abs(diff([num1 num2]))<tol;
   else
      if (nargin<4 || isempty(fh))
         fh = @(x1,x2) abs(x1-x2); 
      end
      if isType(fh, 'string') && any(strcmpi(fh,{'proportion','proport','prop','percent','percentage','perc'}))
         fh  = @(x,y) abs(x./y - 1);
         tol = tol/100;
      end
      
      if ~isscalar(num1) || ~isscalar(num2)
         if isscalar(num1), num1 = ones(size(num2))*num1; end
         if isscalar(num2), num2 = ones(size(num1))*num2; end
         if ~all( size(num1) == size(num2) )
            str = sprintf('\nWhen both inputs are vectors their sizes must match, exiting\n\n');
            cprintf('Keywords*', str);
            return;
         end
         for ii=1:length( toVec(num1) )
            same(ii) = compareFloats( num1(ii), num2(ii), tol, fh );
         end
         return;
      end
      
      % in this case only 1 or neither of the inputs are scalar
      same = fh(num1, num2)<tol | fh(num2, num1)<tol;
   end
      
end