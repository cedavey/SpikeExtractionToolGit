% [cov_x_Yvec,cov_Yvec_x] = cov_x_Yvec(x, Y)
% Calculates covariance between a scalar, x, and a vector Y. The regular 
% cov function requires x and y to have identical sizes if they're
% multidimensional. This function allows x to be a scalar (and therefore
% the sample input is a vector) and Y to be a vector of M dimensions (and
% therefore the sample input is a matrix). Let N be the number of samples,
% which must match for all time series.
% Inputs:
%   x - vector time series of scalar samples; x = [x_1, x_2, ..., x_N]^T
%   Y - matrix containing 2+ time series; Y = [y_1, y_2, ..., y_N]^T, where
%           y_1 = [y_11, y_12, ..., y_1N]^T
%                               :
%           y_1 = [y_M1, y_M2, ..., y_MN]^T, for M dimensions
% Outputs:
%   cov_x_yvec = cov(x, Y), a row vector of size 1xM
%   cov_yvec_x = cov(Y, x), a col vector of size Mx1
function [cov_x_Yvec,cov_Yvec_x] = cov_x_Yvec(x, Y)
    %######## Sanity check of inputs ########%
    if ~isType(x, 'vector') || ~isnumeric(x)
        error('Input error: x must be a vector of numbers'); 
    end
    if ~isType(Y, 'vector') && ~isType(Y, 'matrix')
        error('Input error: Y must be a matrix of numbers'); 
    end
    
    % Make sure that each time series is a col vector
    if size(x,1)  < size(x,2), x = x'; end
    if size(Y,1) ~= size(x,1), Y = Y'; end % num samples must match
    if size(Y,1) ~= size(x,1) % 
        % After transpose if samples still don't match, throw error
        error('Input error: x & Y sample sizes must match');
    end
    
    % Concatenate x & Y, then calculate the cov matrix, wh has form;
    %   |  var(x)  cov(x,Y) |
    %   | cov(Y,x)  var(Y)  |
    cov_matrix = cov([x Y]);
    cov_x_Yvec = cov_matrix(1,2:end);
    cov_Yvec_x = cov_matrix(2:end,1);
end