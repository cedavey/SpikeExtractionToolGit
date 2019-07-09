% EXTRACTSPIKESUSINGKMEANS Classifies identified spikes into different
% clusters and returns each family as a single axon. All the spikes of a
% family are returned in a vector of lenght T (time of the recording).
%
% Syntax: [APspikes, APtimes, APfamily] = extractSpikesUsingKmeans(tseries, method_params, debug);
%
% Created by Artemio - 24/June/2019
function [APspikes, APtimes, APfamilies] = extractSpikesUsingKmeans(tseries, method_params, debugOption)
error('This function hasn''t been finished');
   [spikes, stimes, ~] = getSpikesByThresholding(tseries, method_params);
   [APfamilies, APtimes] = getKmeansClusters(spikes,stimes, debugOption, [10 100]);
   
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
   if ~params.plot_spikes_consecutively.value
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
               fam_tseries(sind,ff) = fam{ff}.spikes(:);
               fam_stimes{ff}  = fam{ff}.stimes(peakN,:);
            end
            APspikes{ti} = fam_tseries;
            APtimes{ti}  = fam_stimes;
            
         catch ME
            % running outta memory so treat each family separately so it's
            % event based rather than disrete time
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
      maxspikes = max( cellfun(@(template) max( cellfun(@(family) length(family.stimes), template ) ), APfamilies ) );
      for ti=1:nAP % for each AP template
         fam = APfamilies{ti}; 
         nF  = length(fam);
         fam_tseries = zeros(maxspikes * nT, nF); % max spikes*spike length x num families
         fam_stimes  = cell(nF, 1);
         for ff=1:nF % for each family within the AP templates
            nspikes = size(fam{ff}.stimes, 2); % number of spikes for this family
            % populate family timeseries by running spikes one after the other
            fam_tseries(1:(nspikes*nT), ff) = toVec(fam{ff}.spikes(:));
            fam_stimes{ff} = fam{ff}.stimes(peakN,:);
         end
         APspikes{ti} = fam_tseries;
         APtimes{ti}  = fam_stimes;
      end
      
   end
   
end