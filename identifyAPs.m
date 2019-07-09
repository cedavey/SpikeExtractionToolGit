% [APs, time] = identifyAPs(tseries, method, method_params)
% Extract AP families from voltage timeseries using thresholding for
% positive and negative peaks, as well as minimum duration positive and
% negative components. Require APs to have user specified similarity to be
% considered part of the same family.
function [APtemplates, componentAPs] = identifyAPs(tseries, method, params, varargin)
   if nargin > 3, debug = varargin{1}; else, debug = 'none'; end
   if nargin==0
      help identifyAPs;
      return;
   end 
   APtemplates = []; componentAPs = [];
   matchtype  = params.match_type.value; % statistic for matching (cov/corr)
   switch lower(method)
      case 'threshold'
         [spikes, stimes, sindices] = getSpikesByThresholding(tseries, params);
      case 'wavelets'
         % now most of the signal is 0's we can use spike extraction tool
         [spikes, stimes] = getWaveletSpikes(tseries, params);
      case 'k means'
         [spikes, stimes, ~] = getSpikesByThresholding(tseries, params);
         [componentAPs, APtemplates] = getKmeansClusters(spikes,stimes, varargin);
         return
      otherwise 
         str = 'No such method for identifying AP templates, ...exiting';
         displayErrorMsg(str);
         return;
   end
   if isempty( spikes ), return; end
   if iscell( stimes ), stimes = cell2mat(stimes(:)); end
   
   [APtemplates, Nsamples, componentAPs] = identifyUniqueAPs(spikes, matchtype, ...
                                                params.match_similarity.value,...
                                                params.normalise_aps.value);
   
   % user may request to ditch AP templates made from only a few spikes
   if params.remove_small_templates.value>0
      remove = Nsamples <= params.remove_small_templates.value;
      APtemplates(:,remove)= [];
      Nsamples(remove)     = [];
      componentAPs(remove) = [];
   end
   

end



% Old method for making spikes the same length but not aligned at peaks
   %% Make all spikes the same duration so can calc cov/corr with potential spikes later
   % get mean spike length 
%    numSamples = cellfun(@length, spikes); 
%    meanLength = double(ceil(mean(numSamples)));
%    outlier    = numSamples < (meanLength - 3*std(numSamples)) | (meanLength + std(numSamples)) < numSamples; 
%    spikes(outlier)     = []; 
%    stimes(outlier)     = [];
%    numSamples(outlier) = [];
%    numSpikes           = length(spikes);
%    
%    % for spikes that are shorter or longer than this, truncate or grab
%    % extra so that all spikes are the same duration
%    tooLong    = numSamples > meanLength; 
%    tooShort   = numSamples < meanLength;
%    
%    % too long --> truncate spike & time vectors
%    spikes(tooLong) = cellfun(@(sp) sp(1:meanLength), spikes(tooLong), 'uni', 0); 
%    stimes(tooLong) = cellfun(@(sp) sp(1:meanLength), stimes(tooLong), 'uni', 0); 
%    
%    % too short --> grab extra, so zero pad each to avoid error msgs
%    data = [tseries.data; zeros(round(meanLength/tseries.dt), 1)];
%    time = [tseries.time; zeros(round(meanLength/tseries.dt), 1)];
%    % funky stuff going on with single datatype so convert
%    sind = @(i) round(double(stimes{i}(1)) / double(tseries.dt));
%    eind = @(i) sind(i) + meanLength - 1;
%    spikes(tooShort) = arrayfun(@(i) data(sind(i):eind(i)), find(tooShort), 'uni', 0);
%    stimes(tooShort) = arrayfun(@(i) time(sind(i):eind(i)), find(tooShort), 'uni', 0);
%    
%    tmp_spikes = reshape( cell2mat(spikes), [meanLength numSpikes]);
