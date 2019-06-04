%  [stat, xn, index] = movingAverageVariance(x, winsize, skip, fn, circ)
% Calculates a moving average of the variance of x using window size w. 
% Where x is a matrix the moving average is calculated for each column. 
% I.e. the format of x must be samples x variables. 
%
% Functions other than variance can by applied by providing a function
% handle in the 'fh' variable. 
%
% Can assume stationarity and calculate a circular moving average. If circ
% is false then the output timeseries will be shorter than the input, with
% length depending on window size. Defaults to false.
%
% Theta allows different weights when calculating the moving average.
%
% Inputs:
%   x     - input data with format samples x dimensions/variable
%   w     - window length for averaging over
%   skip  - number of samples to skip between each window calculation
%   fh    - function handle for summarising data within window if different
%           to variance. 
%   circ  - boolean to use circular moving average or not
%   time  - if input a time vector the output index vector is in time,
%           rather than simply index number
%
% Outputs:
%   stat  - requested statistic calculated across the timeseries (vector)
%   xn    - timeseries divided/normalised by statistic at each sample 
%   index - new timeseries index vector, e.g. [25, 50, 75, 100, ...]
function varargout = movingAverageVariance(x, w, skip, fh, circ, time)
    if nargin==0
        help movingAverageVariance;
        return;
    end
    if nargin<2, w = ceil(length(x)/20); end
    if nargin<3, skip = 1;               end
    if nargin<4, fh   = @var;            end
    if nargin<5, circ = false;           end
    if nargin<6, time = [];              end
    
    % for more complicated functions we need x to be type double
    x = double(x);
    
    szx = size(x);
    changeShape = 0; % haven't changed the shape of x yet

    % if x is a row vector change to col, but change back later
    if isvector(x) && (size(x,1)<size(x,2))
        x = x';  % should use changeShape code instead? Twas added later...
        isRowVec=1;
    else
        isRowVec=0;
        if isDim(x)>2
            x = reshape(x,[prod(szx(1:end-1)) szx(end)])';
            changeShape = 1;
        end
    end
    
    Nin = size(x,1); % number of input samples
    
    % if window size equals number of samples, get the mean
    if isscalar(w) && w>=Nin
        stat = repmat(mean(x,1),[Nin 1]);
        if isRowVec
            stat = stat';
        end
        return
    end
    
    % size of output depends on whether the moving average is circular, &
    % if you skip a chunk of samples each time you move the window forward
    % (if not circular, start taking average after already have a window
    % full of samples)
    if circ
       Nout  = floor(Nin / skip );  % max number we can output, but last window might not fit
    else
       Nout  = floor( (Nin - w) / skip );
       Nout  = Nout - ceil( (w - mod(Nin, w))/skip/2 ); % subtract samples from remaining window that won't fit
    end
    stat  = zeros(Nout, size(x,2));
    if nargout >= 2
       xn = zeros(size(x));
    end
    if nargout >= 3
       Nout = zeros(Nout, size(x,2));
    end
    out = 1; % output index
    if circ
        si= ceil(w/2);        % start index for using full window size
        ei= Nin-floor(w/2)+1; % end index for using full window size
        wl= ceil(w/2);        % get w from left & right of current sample
        wr= w-wl;             % size of window to right of current sample
        for i=1:skip:si
            stat((i-1)/skip+1,:) = fh([x(1:i+wr,:); x(end-wl+i:end,:)]);
            Nout(out) = i; out=out+1; 
        end
        for i=si+1:skip:(ei-1)
            stat((i-1)/skip+1,:) = fh(x(i-wl:i+wr,:));
            Nout(out) = i; out=out+1; 
        end
        for i=ei:skip:Nin
            stat((i-1)/skip+1,:) = fh([x(i-wl:end,:); x(1:wr-(Nin-i),:)]);
            Nout(out) = i; out=out+1; 
        end

    else
        si= floor(w/2);     % start index for using full window size
        ei= Nin-floor(w/2); % end index for using full window size
        wl= floor(w/2);     % get w from left & right of current sample
        wr= ceil(w/2);      % size of window to right or equal to current sample
        for i=si:skip:ei
            st  = ternaryOp((i-wl+1)<1, 1, i-wl+1); % start at index
            ed  = ternaryOp((i+wr  )>Nin, Nin, i+wr  ); % end at index
            ind = max([round((i-wl-1)/skip+1) 1]);

            switch func2str( fh )
               case {'nanstd','std', 'nanmean','mean', 'nanmode','mode'}
                  stat(ind,:) = fh( x(st:ed,:) );
                  if nargout >= 2
                     xn(st:ed,:) = x(st:ed) ./ stat(ind,:);
                  end

               case {'studt','stdGenStudt'}
                  [sig, t]    = stdGenStudt( x(st:ed,:) ); 
                  stat(ind,:) = sig;
                  if nargout >= 2
                     xn(st:ed,:) = normStudt( x, t );
                  end

               case {'varGenStudt'}
                  [sig, t]    = varGenStudt( x(st:ed,:) ); 
                  stat(ind,:) = sig;
                  if nargout >= 2
                     xn(st:ed,:) = normStudt( x, t );
                  end

               case {'nanvar','var'}
                  stat(ind,:) = fh( x(st:ed,:) );
                  if nargout >= 2
                     xn(st:ed,:) = x(st:ed) ./ sqrt( stat(ind,:) );
                  end
               otherwise
                  stat(ind,:) = fh( x(st:ed,:) );
                  xn(st:ed,:) = x(st:ed) ./ stat(ind,:);
            end
            if nargout >= 3
               Nout(out)= i;
               out = out + 1;
            end
        end
    end

   if nargout>=1, varargout{1} = stat; end
   if nargout>=2, varargout{2} = xn;   end
   if nargout>=3
      if isempty( time )
         varargout{3} = Nout; % floor(w/2):skip:(Nin-floor(w/2)); 
      else
         varargout{3} = time( Nout ); 
      end
   end


   if isRowVec
       stat = stat';
   end
   if changeShape
       stat = reshape(stat',szx);
   end
    
end













