%   rho = acf(X, M, [do_plot])
% Matlab function to calculate autocorrelation function (ACF)
% of a time series. The ACF is calculated for various time lags.
% Input:
%	X (vector)       - time series data
%	M (integer)	     - maximum time lag (max=T-1)
%	do_plot (binary) - 1 if you want the ACF plotted, 0 o/wise
% Output:
%	rho (vector)     - gives ACF for each time lag between 1 & M
function [rho] = acf(x, varargin)
    % deal with variable number of inputs - set default values, then write
    % over with user initialised variables
    optargs = {10 0}; % optargs = {M(max time lag) do_plot}
    numvarargs = length(varargin);
    if numvarargs>2
        error('acf:TooManyInputs','requires at most 2 optional args');
    end
    optargs(1:numvarargs) = varargin;
    [M, do_plot] = optargs{:};
    
    if isempty(M)
        M=length(x)-1;
    elseif M>=length(x)
        error('Max lag is T-1 (T=num time pts)');
    end
    x  = x(:); % force x to column vector
    x  = x - mean(x);
    T  = length(x); % number of samples
    rho=conv(flipud(x),x);
    rho=rho(T:end);
    den=[T:-1:1]';
    rho=rho./den;
	if rho(1)==0
		rho(:) = 0;
	else
	    rho=rho/rho(1);
	end

    if isscalar(M) && M~=0
        rho = rho(1:M+1); % M+1 cuz' 1st elmt is lag 0
        M = 0:M;
    elseif isvector(M)
        rho = rho(M+1);
    else % M = 0 (just included for completeness)
        rho = rho(1); 
    end
    
    if do_plot
       figure
       % don't plot rho(1) since it's always 1 (& maybe >> so lose resolution)
%         bar(1:M, rho(2:end)); 
        bar(M, rho); % plotting rho(1) for article consistency
        xlim([min(M) max(M)]);
    end % end do plot

end % end function
