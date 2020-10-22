% [mavdata, mavtime] = movingAverageUnevenSamples( time, data, avtime, avfn )
%
% Calculate moving average of data timeseries where time samples are not
% evenly spaced. Returns averaged data timeseries & the time points they
% occur
%
% Inputs:
%  time    - unevenly spaced time vector
%  data    - data at time points given in time vector
%  avtime  - duration of averaging period, in same units as the time vector
%  skiptime- after each avging jump ahead this amount, in same units as time
%  avfn    - function used for averaging (defaults to @mean)
% Outputs:
%  mavdata - moving average data vector
%  mavtime - time points for moving average data vector (assumed to be
%            causal so if averaging over a window, time is the endpoint of
%            the window)
function [mavdata, mavtime] = movingAverageUnevenSamples( time, data, avtime, skiptime, avfn )
   % TO DO: use matrix operations, take out the loop, to improve efficiency
   if nargin<4 || isempty( skiptime ), skiptime = avtime; end
   if nargin<5 || isempty( avfn ), avfn = @mean; end
   mavdata  = zeros( round( (time(end) - time(1)) / skiptime ) - 1, 1 ); % -1 cuz start at avtime, not skiptime
   mavtime  = zeros( size( mavdata ) );
   
   currtime = time(1) + avtime;
   avgnum   = 1;
   while currtime < time(end)
      prevtime = currtime - avtime; 
      tind     = prevtime <= time & time <= currtime; 
      mavdata(avgnum) = avfn( data( tind ) );
      mavtime(avgnum) = currtime; 
      avgnum   = avgnum + 1;
      
      currtime = currtime + skiptime;
   end
end