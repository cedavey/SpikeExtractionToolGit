% EXTRACTSPIKESUSINGKMEANS Classifies identified spikes into different
% clusters and returns each family as a single axon. All the spikes of a
% family are returned in a vector of lenght T (time of the recording).
%
% Syntax: [APspikes, APtimes, APfamily] = extractSpikesUsingKmeans(tseries, method_params, debug);
%
% Created by Artemio - 24/June/2019
function [APspikes, APtimes, APfamilies] = extractSpikesUsingKmeans(tseries, method_params, debugOption)
   [spikes, stimes, ~] = getSpikesByThresholding(tseries, method_params);
   [APfamilies, APtimes] = getKmeansClusters(spikes,stimes, debugOption, [1 20]);
   
   % Get peakN (index of peak value
   [~, idx] = max(APfamilies{1}{1}.spikes');
   peakN = round(mean(idx));
   % get rid of empty families, just in case they snuck in there!
   nf = cellfun( @length, APfamilies );
   loner = nf <= 0; % get empty families
   APfamilies( loner ) = [];
   APtimes ( loner ) = [];
   if isempty( APfamilies )
      displayErrorMsg( 'No axon families found with this parameter configuration' );
      return; 
   end
   nAP = length( APfamilies ); % update num families after removing loners
   
   APspikes = cell(1,nAP);
   % Each valid spike has been allocated to an APfamiliy. Now for each 
   % AP family we need to construct a timeseries
   if ~method_params.plot_spikes_consecutively.value
      for ti=1:nAP % for each AP template
         fam = APfamilies{ti}; 
         nF  = length(fam);
         
         try
            % try allocating a full tseries for each family in template
         	fam_tseries = zeros(length(tseries.data),nF);
         	fam_stimes  = cell(nF, 1);
            
            for ff=1:nF % for each family within the AP templates
               % if we're creating the timeseries extract spikes at
               % appropriate times & copy into voltage timeseries
               sind = round(fam{ff}.stimes(:) / tseries.dt);
               sind(sind == 0) = [];
               fam_tseries(abs(sind),ff) = fam{ff}.spikes(:);
               if size(fam{ff}.stimes,1) == 1
                  fam_stimes{ff}  = fam{ff}.stimes(1,peakN);
               else
%                   fam_stimes{ff}  = fam{ff}.stimes(peakN,:);
                  fam_stimes{ff}  = fam{ff}.stimes(peakN,1);
               end
            end
            APspikes{ti} = fam_tseries;
            APtimes{ti}  = fam_stimes;
            
         catch ME
            % running outta memory so treat each family separately so it's
            % event based rather than discrete time
            for ff=1:nF 
               % for each family within the AP templates, extract the times of
               % spike peaks
               fam_stimes{ff}  = fam{ff}.stimes(peakN,:);
            end
            fam{1}.time     = tseries.time;
            APspikes{ti} = fam;
            APtimes{ti}  = fam_stimes;
         end
      end
      
   % squish spikes together so they can be seen as one long series of spikes
   else
      % determine which family has the most spikes, because we'll make all
      % timeseries accommodate that length (currently time limits are
      % common to whole dataset)
      maxspikes = 0;
      for i = 1:size(APfamilies)
         if size(APfamilies{i}{1}.spikes,1) > maxspikes
            maxspikes = size(APfamilies{i}{1}.spikes,1);
         end
      end
      spikesize = maxspikes * size(APfamilies{1}{1}.spikes,2);
      APspikes = cell(size(APtimes'));
      APtimes = cell(size(APspikes));
      for ti=1:nAP % for each AP template
         fam = APfamilies{ti}{1}; 
         APspikes{ti} = zeros(spikesize, 1);
         fam_tseries = toVec(fam.spikes');
         APspikes{ti}(1: length(fam_tseries)) = fam_tseries;
         APtimes{1,ti} = cell(1);
         APtimes{1,ti}  = {fam.stimes(:,peakN)};
      end
      
   end
   
end