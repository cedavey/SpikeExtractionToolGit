function [spikes, stimes] = getWaveletSpikes(tseries, params)
   spikes     = []; stimes = [];
   [tseries, err] = filterWithWavelet(tseries, params);
   voltage    = tseries.data; time = tseries.time; 
   nT         = length( voltage );
   dt         = tseries.dt;
   minsp      = 2e-3; % min spike length in seconds
   minlength  = ceil( minsp / dt);
   clear tseries;
   if err, return; end
   
   % go through wavelet times series and extract spikes, then we'll
   % identify unique templates within the spikes 
   
   % to do: need to find a way of identifying multiple spikes in a row, cuz
   % at the mo it gives 1 v.long spike, wh is ditched for being too long
   % (perhaps keep wavelet instead of substituting voltage back in)
   isspike    = voltage ~= 0;
   % to ensure we identify a spike at v.begining & end, prepend & append 0s
   spikestart = find( diff([0; isspike; 0])== 1 );
   spikeend   = find( diff([0; isspike; 0])==-1 );  
   splength   = spikeend - spikestart + 1;
   % ditch v.long or v.short spikes (get rid of long 1st cuz they skew mean)
   ditch      = (mean(splength) + 3*std(splength)) < splength;
   spikestart = spikestart( ~ditch );
   spikeend   = spikeend( ~ditch );
   % short spikes are only 1 std dev from mean cuz it's heavy tailed at
   % long end so if minus 3 or 2 times std dev it becomes -ve
   splength   = spikeend - spikestart + 1;
   ditch      = splength < minlength;
   spikestart = spikestart( ~ditch );
   spikeend   = spikeend( ~ditch );
   splength   = spikeend - spikestart + 1;
   
   % go through long spikes and make sure it's not consecutive valid spikes
   % (removed cuz the spikes seemed too short, wh is perhaps why they were
   % consecutive in the first place?)
%    for di=1:ditch
%       tmpvolt = voltage( spikestart(ditch(di)) : spikeend(ditch(di)) );
%       dvolt   = diff( tmpvolt );
%       posv    = tmpvolt>0;
%       startpos= find(diff([0; posv; 0])==1); % location of change from 0 to 1
%       startneg= find(diff([0; posv; 0])==-1);
%       durpos  = startneg - startpos;
%    end
   % make sure last spike doesn't go beyond voltage length cuz of appended 0s
   if spikeend(end) > length(voltage), spikeend(end) = length(voltage); end
   
   spikes     = arrayfun(@(i) voltage( spikestart(i):spikeend(i) ), ...
                         1:length(spikestart), 'un', 0);
   stimes     = arrayfun(@(i) toVec(time(spikestart(i):(spikeend(i)))), ...
                              1:length(spikestart), 'uniformoutput', false); 
   pospeak    = cellfun(@(sp) find(sp==max(sp),1,'first'), spikes);
   % also get rid of spikes where peak -ve happens before peak +ve
   negafter   = min(6, minlength);
   negpeak    = cellfun(@(sp) find(sp(negafter:end)==min(sp(negafter:end)),1,'first')+negafter-1, spikes);
   ditch      = pospeak > negpeak;
   spikes     = spikes( ~ditch );
   stimes     = stimes( ~ditch );
   pospeak    = cellfun(@(sp) find(sp==max(sp),1,'first'), spikes);
   
   tbefore    = ceil(mean(pospeak));
   splength   = cellfun(@length, spikes);
   splength   = ceil(mean(splength));
   tafter     = ceil(mean(splength) - tbefore);
   % get mean time before & after spike peak - extend voltage so don't go over the edge
   voltage    = [voltage; zeros(tafter,1)];
   % if first spike peaks before avg pre-peak time, add zeros to it
   if round(stimes{1}(pospeak(1))/dt) < tbefore
      voltage = [zeros(tbefore,1); voltage];
      shift   = tbefore;
   else
      shift   = 0;
   end
   sp         = arrayfun(@(i) voltage(round(stimes{i}(pospeak(i))/dt)+shift-tbefore+1 : round(stimes{i}(pospeak(i))/dt)+shift+tafter), ...
                              1:length(spikes),'uniformoutput', false);
   spikes     = cell2mat(sp);
   stimes     = arrayfun(@(i) time(round(stimes{i}(pospeak(i))/dt)+shift-tbefore+1 : round(stimes{i}(pospeak(i))/dt)+shift+tafter), ...
                              1:length(spikes),'uniformoutput', false);
end

function [tseries, err] = filterWithWavelet(tseries, params)
   dt      = tseries.dt;   
   fs      = 1/dt;
   f_min   = params.lower_frequency.value;
   f_max   = params.upper_frequency.value;
   thresh  = params.wavelet_threshold.value;
   % see exploreWavelets.m for explanation of the following
   No      = 8; % number of octaves
   Nv      = 8; % number of voices
   io      = (1:No)'; iv = (1:Nv)'; 
   sfn     = @(io, iv) 2.^(io + iv/Nv);
   scales  = toVec( sfn( repmat(io',[Nv 1]), repmat(iv,[1 No]) ) );
   fscales = scal2frq(scales, 'morl', dt);
   durAP   = 6e-3; % spike typically takes this long 
   if f_max==0, f_max = 1 / durAP / 2; end
   if f_min==0, f_min = 1 / durAP * 2; end
   
   % try running continuous time wavelet transform
   err = true;
   while err && No>=3 && Nv>=4
      try
         % wavelet params: symmetry, & time-bandwidth in number of dt's
         [coeff,f] = cwt(double( tseries.data ), 'morse', fs, ...
                         'NumOctaves', No, 'VoicesPerOctave', Nv, ...
                         'WaveletParameters', [3 6]); 
         if length(unique(coeff(:)))==1
            displayErrorMsg('Wavelet coeffs all zero');
            return;
         end
         err = false;
      catch ME
         % if we fail it's most likely a memory issue
         No      = ceil(No*2/3);  % number of octaves
         Nv      = ceil(Nv*2/3);  % number of voices
         io      = (1:No)'; iv = (1:Nv)'; 
         sfn     = @(io, iv) 2.^(io + iv/Nv);
         scales  = toVec( sfn( repmat(io',[Nv 1]), repmat(iv,[1 No]) ) );
         fscales = scal2frq(scales, 'morl', dt);       
      end
   end
   % if we're here & err is still true, we weren't successful but rather
   % quit trying because number of octaves & voices got too small
   if err
      displayErrorMsg('Memory error, had to give it up sorry');
      return;
   end
   % coeff format is freq/scale x time, where coeff(1,:) is highest freq &
   fmind   = find( (fscales - f_min)<0, 1, 'first');
   fMind   = find( (fscales - f_max)<0, 1, 'first');
   % coeff format is freq x time, with freqs in descending order
   if isempty(fmind), fmind = size(coeff,1); end
   if isempty(fMind), fMind = 1;             end
   
   % icwt will bandpass for us, so just threshold coefficients for now
   % coeffs for each freq are pretty gauss
   thresh = thresh*std(coeff(:));
   coeff( abs(coeff) < thresh ) = 0;

   voltage = icwt( coeff, 'morse', f, [fscales(fmind) fscales(fMind)], 'WaveletParameters', [3 6] );
   if length(unique(voltage))==1
      displayErrorMsg('Voltage signal built from wavelets is empty');
      return;
   end
   % replace voltage with orig values since reconstructed APs may not be as
   % sharp or accurate after removing freqs
   voltage( voltage~=0 ) = tseries.data( voltage~=0 );
   
   tseries.data = voltage(:);
end