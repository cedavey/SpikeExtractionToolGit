% true = isDim(X, dim) 
% isDim returns true if the input X is the asserted size, dim.
% Returns false otherwise. If dim is not provided the dimension of X is
% returned.
% Inputs:
%    X    - is the X to determine the size of
%    dim  - asserted size of X (can be vector - returns boolean for each value)
% Outputs: true = isDim(X, dim) 
%   true  - returns 1 if X is of size dim, 0 otherwise
%   X_    - same as X but with singleton dimensions removed using squeeze
%
% Outputs: dim = isDim(X)
%   dim   - returns the dimension of X
%   X_    - same as X but with singleton dimensions removed using squeeze
function [ok, X] = isDim(X, dim)
%     X = squeeze(X);
    sz = length(size(X));

    if nargin==2
        ok = 0;
        % test for vectors separately because vectors have 2 dims w above test
        if isvector(X)
            if dim==1
                ok = 1;
            end
            return
        end
        ok = (any(sz==dim)); % dim may be a vector
    else
      % test for vectors separately because vectors have 2 dims w above test
      if isvector(X), ok = 1; return; end;
      ok = sz;
    end
end
