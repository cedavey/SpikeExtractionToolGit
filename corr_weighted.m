function [C] = corr_weighted(W,X,Y)
    % corr_weighted Description: Performs a weighted correlation between
    % each column of X, or X and Y using weights W, where W is the same
    % size as X and Y
    % ______________
    %[C] = corr_weighted(W,X,Y) : correlate all elements of X with all
    %   elements of Y, with weighting W
    %[C] = corr_weighted(W,X) pairwise correlation of columns of x with
    %weighting W
    
    % Inputs:
    %   - W: weight matrix for each index comparison, same size as X
    %   - X: matrix to correlate. if Y is input, then we use each element
    %   weighted by W, if only X, then we correlate each column with each
    %   other column. In the X-Y correlation case, we have one set of
    %   weights that corresponds to the elements of each set, in X
    %   columnwise case we use the weights corresponding to each column to
    %   weight the columns separately
    %   - Y:
    % Outputs:
    %   - C:
    % ______________
    % Author: Kathrine Clarke 
    % Date: 29-Sep-2023
    
    if nargin == 3
        C = Wcorr(W,W,X,Y);
    else
        C = zeros(size(X,2));
        for col1 = 1:size(X,2)
           for col2 = 1:size(X,2)
                C(col1,col2) = Wcorr(W(:,col1),W(:,col2),X(:,col1),X(:,col2));
           end
        end
    end

end
 
function C = Wcorr(Wx,Wy,X,Y)
    Wxy = max(Wx,Wy);
    Wxy_sum = sum(Wxy,'all');
    mX = sum(Wxy.*X,'all')./Wxy_sum;
    mY = sum(Wxy.*Y,'all')./Wxy_sum;
    sX = sum(Wxy.*(X-mX).^2,'all')./Wxy_sum;
    sY = sum(Wxy.*(Y-mY).^2,'all')./Wxy_sum;
    sXY = sum(Wxy.*(X-mX).*(Y-mY),'all')./Wxy_sum;
    C = sXY./sqrt(sX.*sY);
end