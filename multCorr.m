% [R, R2 Rprob] = multCorr(x, Y)
% Calculates multiple correlation between a scalar, x, and a vector Y.
%            R^2 = (cov(x,Y).var(Y)^-1.cov(Y,x))/var(x)
% Inputs (assuming N samples of each variable):
%   x    - time series (column vector) containing samples of a scalar, x_t
%          size: N x 1
%   Y    - matrix of 2+ time series, [y_1t,y_2t,...,y_Mt] at each time point
%          size: N x M, where N is num samples, & M is num predictors
% Outputs:
%   R    - multiple correlation coefficient (vector version of correlation)
%  R2    - coefficient of multiple determination. R2 = R*R
%  Rprob - probability of R compared to the null hypothesis (no predictors)
function [R, R2, Rprob] = multCorr(x, Y)
    %######## Sanity check of inputs ########%
    if ~isType(x, 'vector')
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

    % we only want columns of Y that contain unique values - if it's a
    % column of constant values for obtaining the mean it'll give a 0
    % covariance & therefore a singular covariance matrix
    j=1;
    for i=1:size(Y,2), 
        if length(unique(Y(:,j)))==1
            Y(:,j)=[];
            j=j-1;
            if isempty(Y)
                R=NaN; R2=NaN; Rprob=NaN;
                return;
            end
        end
        j=j+1;
    end
    %##### R^2 = (cov(x,Y).var(Y)^-1.cov(Y,x))/var(x) #####%
    [cov_xY,cov_Yx] = cov_x_Yvec(x, Y); % cov(x,Y) & cov(Y,x)
    warning off
    varY_inv = cov(Y)^-1;               % var(Y)^-1
    if any(isinf(varY_inv(:)))
        varY_inv = pinv(cov(Y));
    end
    warning on
    varx = var(x);                      % var(x)
    R2 = (cov_xY*varY_inv*cov_Yx)/varx; % R^2
    R = sqrt(R2);                       % R = sqrt(R^2)
    m = size(Y,2); % number of predictors
    M = size(x,1); % number of samples
    Rtest = R2/(1-R2)*M/m;
    Rprob = fcdf(Rtest,m,M);
end