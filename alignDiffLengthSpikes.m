% [ sp1, sp2, time1, time2 ] = alignDiffLengthSpikes( sp1, sp2, time1, time2, varargin )
%
% For spikes that are different lengths, this function aligns them at the
% peak and then gets the min/max number of samples of either spikes before the 
% peak, and the min/max number of samples of either spike after the peak, to 
% ensure both spikes have valid samples when aligned. If time vectors are
% provided they are adjusted to match the size of the resulting spikes. If
% an input time is from a family (because one of the input spikes is a
% template), it will be a matrix representing spike times for all spikes in
% the template, and all spikes times will be adjusted accordingly. 
%
% varargin allows other spikes to be added that will not be used to
% determine where the peaks match, but will be adjusted in line with the
% spike 1, for example where sp1 is a mean template, and the varargin 
% spikes are spikes in the template family (only works for sp1). 
function [ sp1, sp2, time1, time2, varargout ] = alignDiffLengthSpikes( sp1, sp2, time1, time2, shorten, varargin )
   markerfn = @getMaxInd; % @getZeroCrossingInd; %
   if nargin<5 || isempty(shorten), shorten = false; end
   
   % if spikes & templates can have diff lengths, use the smaller length
   % to calculate the match, but centre both around the peak 
   varargout = cell(0);
   if nargin>2 && ~isempty(time1) && ~isempty(time2)
      havetime = true; 
   else
      havetime = false; 
      time1    = [];
      time2    = [];
   end

   % remove any NaNs to begin with (sometimes added to align)
   sp1_nan  = isnan(sp1);
   sp2_nan  = isnan(sp2);
   sp1(sp1_nan) = 0;
   sp2(sp2_nan) = 0;
   if havetime
      dt    = mean(diff(time1));
      time1(sp1_nan) = 0;
      time2(sp2_nan) = 0;
   end
   nSp1     = size( sp1, 1 );
   nSp2     = size( sp2, 1 ); % length of each template
   peakSp1  = round( mean( markerfn( sp1 ) ) ); % in case of many spikes
   peakSp2  = round( mean( markerfn( sp2 ) ) );

   if shorten
      tbefore = min( peakSp2, peakSp1 );
      tafter  = min( nSp2 - peakSp2, nSp1 - peakSp1 );
   else
      tbefore = max( peakSp2, peakSp1 );
      tafter  = max( nSp2 - peakSp2, nSp1 - peakSp1 );
   end

   [sp1, time1] = adjustSpike( shorten, sp1, tbefore, tafter, peakSp1, time1 );
   [sp2, time2] = adjustSpike( shorten, sp2, tbefore, tafter, peakSp2, time2 );
   while ~isempty(varargin)
      varargout{end+1} = adjustSpike( shorten, varargin{:}, tbefore, tafter, peakSp1 );
      varargin(1) = [];
   end
   

end

% Takes an input spike and its peak (or zero crossing), and the time
% required both before and after the peak. Appends 0s if the spike has
% insufficient samples before or after the peak (or zero crossing).
% Note: spike can't be too long because the longest spike sets the lengths.
% If time vector provided, expand the time vector for any added zero samples.
% If any additional inputs provided, assume they need to be modified in the
% same way as the input spike

function [sp_resz, time_resz, varargout] = adjustSpike( shorten, sp, tbefore, tafter, peakSp, sptime, varargin )
   % if time before spike peak (or zero crossing) is too short, append 0s
   % (time before includes the peak)

   % if modify spike length, also modify time vector
   if nargin>5 && ~isempty(sptime)
      havetime  = true; 
      time_resz = sptime; % initialise spike time vector
      dt        = sptime(2) - sptime(1);
   else
      havetime  = false;
      time_resz = [];
      dt        = 0;
   end

   % ensure input spike is a column vector, if it is a vector at all, so
   % when modifying the spike we're working in the first dimension
   if isvector(sp)
      iscol = isColVec(sp);
      if ~iscol
         sp = sp'; 
         if havetime, sptime = sptime(:); end
      end
   else
      iscol = false;
   end
       
   N         = size(sp,1);
   sp_resz   = sp;
   varargout = cell(0);
   
   % functions for prepending and appending zeros
   if shorten
      prependTest = @(peak)      tbefore < peak;
      appendTest  = @(peak)      tafter  < (N - peak);
      prependfn   = @(sp, n)     sp(n+1:end,:);
      appendfn    = @(sp, n)     sp(1:end-n,:);
      prependTime = @(time,n,dt) time(n+1:end,:);
      appendTime  = @(time,n,dt) time(1:end-n,:);

   else
      prependTest = @(peak)      tbefore > peakSp;
      appendTest  = @(peak)      tafter  > (N - peak);
      prependfn   = @(sp, n)     [zeros(n,size(sp,2)); sp];
      appendfn    = @(sp, n)     [sp; zeros(n,size(sp,2))];
      prependTime = @(time,n,dt) [ ( bsxfun(@plus, (time(1,:)-(n+1)*dt)', dt*(1:n)) )'; time ];
      appendTime  = @(time,n,dt) [ time; ( bsxfun(@plus, time(end,:)', dt*(1:n)) )'];
   end


   % Prepend before peak (or zero crossing) if its too short, and append
   % after peak (or zero crossing) if its too short. Make the same
   % adjustments to spike time (but by adding time indices), and to any
   % other spikes provided in varargin (e.g., if modify a family template, 
   % we need to modify the family's spikes also
   if prependTest(peakSp)
      prepend   = abs(tbefore - peakSp);   % needs to work for shortening & lengthening
      sp_resz   = prependfn(sp_resz, prepend);
      if shorten
         N      = N - prepend;
         peakSp = peakSp - prepend;
      else
         N      = N + prepend;
         peakSp = peakSp + prepend;
      end
      if havetime
         time_resz = prependTime(time_resz,prepend,dt);
      end
      if ~isempty(varargin)
         while ~isempty(varargin)
            varargout{1} = prependfn( varargin{1}, prepend );
            varargin(1)  = [];
         end
      end
   end
   if appendTest(peakSp)
      append    = abs(tafter - (N-peakSp));   % needs to work for shortening & lengthening
      sp_resz   = appendfn(sp_resz, append);
      if havetime
         time_resz = appendTime(time_resz,append,dt);
      end
      if ~isempty(varargout)
         argout_copy = varargout; % removing inputs to know they're dealt with, so copy first
         while ~isempty(argout_copy)
            varargout{1}   = appendfn( varargout{1}, append );
            argout_copy(1) = [];
         end
      end

   end

   if isvector(sp)
      % assumes sp and sptime are either both column, or both row, vectors
      if ~iscol
         sp_resz = sp_resz';
         if havetime, time_resz = time_resz'; end
      end
   end
end % end lengthening spike







