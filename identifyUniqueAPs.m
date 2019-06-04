% uniqueAPs = identifyUniqueAPs(APs, matchtype, matchthresh)
% 
% Identify a unique set of APs from the action potentials extracted from a
% voltage time series. Calc corr of AP with current set of unique APs &, if
% corr is less than the threshold for all then create a new unique AP.
% Otherwise merge the AP with the unique AP that it has the highest
% correlation with. The order in which the APs are presented to be
% classified is randomised.
%
% Inputs:
%  APs        - all spike action potentials extracted using extractActionPotentials
%              (format: is time x N)
%  corrthresh - correlation of AP with unique set must be less than this to be unique
%              (format: scalar)
%
% Outputs:
%  uniqueAPs - copy of all unique action potential waveforms 
%              (format: time x k, where k is the number of unique APs)
%  Nsamples  - number of APs used in construction of mean AP waveform (kx1)
function [APtemplates, Nsamples, componentAPs] = identifyUniqueAPs(spikes, matchtype, matchthresh, normAPs)
%% TO DO:
% - if matching on corr --> normalise since we don't care about amplitude
   if ~exist('normAPs','var') || isempty(normAPs), normAPs = false; end
   componentAPs = cell(0);
   
   if nargin<2, matchthresh = 0.9; end
   N = size(spikes, 2);
   
   order     = randperm(N); % randomise the order in which the APs are classified
   
   % start AP template families off by putting first spike in 
   ap = spikes(:, order(1));
   if normAPs
      ap = ap/std(ap);
   end
   APtemplates     = ap; 
   componentAPs{1} = ap; % record all constituents of each AP family
   Nsamples        = 1;  % number of AP samples for each waveform

   for i=2:N
      ap = spikes(:,i); 
      
      switch matchtype
         case 'corr'
            rho = corr(APtemplates, ap);
            
         case 'cov'
            % to compare using covariance means that the APs have to have
            % the same amplitude to be considered part of the same family.
            % In order to have a match threshold independent of amplitude,
            % scale by the size of the template's signal. Since cov is
            % voltage squared, scale by variance rather than std dev
            rho = cov_x_Yvec(ap, APtemplates) ./ var(APtemplates); 
      end
      
      % Current AP not correlated with any waveform, so create new unique AP
      if ~any(rho > matchthresh)
         % make a new AP family - mean AP in uniqueAPs, and add to
         % constituent APs of the family in componentAPs, and record number
         % of APs in this unique AP family (i.e. 1 since it just began)
         if normAPs
            ap = ap / std(ap); 
         end
         APtemplates(:,end+1) = ap;
         componentAPs(end+1)  = cell(1);
         componentAPs{end}    = ap;
         Nsamples(end+1)      = 1;
         
      % Get waveform the AP is maximally correlated with, add it to the
      % group & calculate the new mean. If matching on corr we don't care
      % about amplitude so normalise, else big APs will drown out small APs.
      % If matching on cov we do care about AP size so leave as is.
      else
         [~,k] = max(rho);
         if strcmpi( matchtype, 'corr' ) && normAPs
            ap = ap / std(ap); 
         end
         componentAPs{k}(:, end+1) = ap;
               
         APtemplates(:,k) = (APtemplates(:,k)*Nsamples(k) + ap) / (Nsamples(k)+1);

         Nsamples(k) = Nsamples(k) + 1;
      end     
   end
   
   % now we've finished processing spikes, merge APtemplates that are
   % correlated more than the required threshold - AP templates may have
   % become more similar over time, as spikes have been added
   % obviously we change the correlation every time we merge APs, so
   % randomise the order we merge in, & only merge once
   rho        = corr(APtemplates);
   [nT,nAPs]  = size(APtemplates);
   order      = randperm(nAPs, nAPs);
   for i=1:nAPs
      if i>nAPs, break; end % reducing number of APs each time through so check
      jj      = order(i);
      ind     = 1:nAPs;
      rho     = corr( APtemplates(:,jj), APtemplates );
      rho(jj) = 0;
      merge   = rho > matchthresh;
      if any(merge)
         merge(jj) = true; % set to 1 so it's included in the merge
         % when calculating merged AP template, weight avg by num samples
         N = Nsamples(merge); % weight for each template
         APtemplates(:, jj)    = sum(APtemplates(:, merge) .* repmat(N(:)', [nT 1]) / sum(N), 2); 
         componentAPs{jj}      = cell2mat(componentAPs(merge)); 
         Nsamples(jj)          = sum(Nsamples(merge)); % add samples from all matching templates
         
         % now set primary template back to 0 so we don't ditch it
         merge(jj)             = false; 
         APtemplates(:, merge) = [];
         Nsamples(merge)       = [];
         componentAPs(merge)   = [];
         % remove indices of merged AP templates 
         mergeind = find(merge);
         order    = arrayfun(@(ii) ternaryOp(any(ii==mergeind), 0, ii), order);
         order(order==0) = [];
          % udpate order indices by counting how many of the merge
          % templates were sitting before each template
         order    = arrayfun(@(ii) ii - sum(ii > mergeind), order);
         nAPs     = nAPs - sum(merge);
      end
   end
   
   % now order AP families by how many APs went into them
   [~, order]   = sort(Nsamples, 'descend');
   APtemplates  = APtemplates(:, order);
   componentAPs = componentAPs(order);
   Nsamples     = Nsamples(order);
end

