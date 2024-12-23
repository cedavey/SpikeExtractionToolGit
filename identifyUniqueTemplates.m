% [templates,alignment_ind_templates, Nsamples, templateAPs]  = identifyUniqueTemplates(spikes,alignment_inds, matchtype,matchthresh)
% 
% Identify a unique set of AP templates from the spikes extracted from a
% voltage time series. Calc corr of AP with current set of AP templates &, 
% if corr is less than the threshold for all, indicating the the spike 
% doesn't match any of the templates very well, then create a new AP 
% template. Otherwise, if the spike matches one or more templates sufficiently 
% well, merge the spike with the AP template that it has the highest
% correlation with. If the spike is well correlated with multiple templates, 
% then other features are used to decide, such as the similarity in spike 
% and template peaks or minima. 
% 
% The order in which the APs are presented to be classified is randomised, 
% to avoid the order in which spikes are presented biasing the results.
% This is necessary because when spikes are merged with a template, the
% template shape changes proportionally to how many spikes have contributed
% to creating the template. 
%
% Note that if spikes is a matrix of time x spike_num, then all spikes
% are treated as having the same length. If spikes is a cell array, then
% spikes are not treated as having a common length. Where the spikes
% and template lengths are different, correlation is calculated using the
% shorter length, but merging uses the entire spike, where an ensemble mean
% is taken across each sample within the spike using only the spikes with a
% value for each sample
%
% Inputs:
%  spikes       - all spike action potentials extracted using extractActionPotentials
%                (format: is time x num_spikes, OR cell array with one
%                cell for each spike, and cell format voltage x 1
%  alignment_inds - vector of alignment inds same size as spikes (for
%  variable length code only, leave empty if using fixed length
%  matchtype    - 'corr' or 'cov': use corr to match only on shape, cov to 
%                 match on shape & size (format: string, defaults" 'corr')
%  matchthresh  - correlation of spike with template must be less than this 
%                 to be considered unique & a new template initiated
%                (format: float scalar, default: 0.95)
%  normAPs      - if true normalise AP templates for future spike matching
%                (format: boolean, default: true)
%
% Outputs:
%  templates    - copy of all unique action potential waveforms 
%                 (format: time x k, where k is the number of unique APs      )
%              OR (depending on whether the input 'spikes' is a cell or matrix)
%                 (format: cell x k, and time x 1 within each cell, but the   )
%                  duration of time changes for each template                 )

%  alignment_ind_templates - the alignment index for each template (vector)
%  Nsamples     - number of APs used in construction of mean AP waveform (kx1)
%                 (cell array of matrices, 1 cell per template, and a
%                 single vector of the template, but templates have diff
%                 length vectors)
%  componentAPs - keeps a list of all spikes contributing to the templates 
%                 AP (cell array of matrices, 1 cell per template, and 
%                 matrix size template_length x num_spikes)
function [templates,alignment_ind_templates, Nsamples, templateAPs] = identifyUniqueTemplates(spikes,alignment_inds, varargin)
%% TO DO:
% - if matching on corr --> normalise since we don't care about amplitude
   optargs = {'cov', '0.9', []}; 
   nargs   = length( varargin );
   optargs(1:nargs) = varargin(:);
   [matchtype, matchthresh, normAPs] = optargs{:};
   if isempty( normAPs ), normAPs = ternaryOp( matchtype=='cov', false, true ); end
   
   % if spikes are in cells then each spike has a different length, else if
   % in a matrix then they've been forced to a common length, calculated
   % from the mean spike length of all spikes
   if iscell( spikes ), uniqueAPLength = true; else, uniqueAPLength = false; end
   N = size(spikes, 2);   
   order = randperm(N); % randomise the order in which the APs are classified
   
   % start AP template families off by putting first spike in a template
   if uniqueAPLength
      sp = spikes{ order(1) };
   else
      sp = spikes(:, order(1));
   end
   if normAPs
      sp = sp/std(sp);
   end
   % templates are unique in length --> use cell, common length --> use
   % matrix
   templates      = ternaryOp( uniqueAPLength, {sp}, sp );
   [~, tempPeak ] = max( sp );
   templateAPs{1} = sp; % record all constituents of each AP family
   alignment_ind_templates(1) = alignment_inds(order(1)); %start template alignment indexes with the spike alignment index
   Nsamples       = 1;  % number of AP samples for each waveform
   Ntemplates     = 1;  % record the number of templates we have

   % get match function
   switch matchtype
      case 'corr'
         % KAT: this doesnt work when we don't use cells, because xcorr
         % doesnt take multiple columns - need to change this function
         matchfn = @( sp, temp ) corr( temp, sp );
         matchfn = @( sp, temp ) max(xcorr(temp, sp, 'normalized', 5)); % consider max 5 lags

      case 'cov'
         % to compare using covariance means that the APs have to have
         % the same amplitude to be considered part of the same family.
         % In order to have a match threshold independent of amplitude,
         % scale by the size of the template's signal. Since cov is
         % voltage squared, scale by variance rather than std dev
         matchfn = @( sp, temp ) cov_x_Yvec( sp, temp ) ./ var( temp );

   end
   % make a fn to extract spike based on matrix or cell array
   getspike = ternaryOp( uniqueAPLength, @(i) spikes{i}, @(i) spikes(:,i) );
   
   cprintf( 'keywords', '\tIdentifying unique family templates\n');
   for ii=2:N
      i = order(ii);
      
      sp  = getspike(i); 
      rho = spikeCorrWithTemplates( sp, templates, alignment_inds(i),alignment_ind_templates, matchfn );
      
      % if rho is empty then the spike wasn't able to be matched with
      % templates because less than half of the spike was overlapping in
      % pre & post peak time periods
      if ~isempty( rho )
         % Current spike not correlated with any template, so create new one
         if ~any(rho > matchthresh)
            % make a new AP family - mean AP in uniqueAPs, and add to
            % constituent APs of the family in componentAPs, and record number
            % of APs in this unique AP family (i.e. 1 since it just began)
            if normAPs
               sp = sp / std(sp); 
            end
            Ntemplates  = Ntemplates + 1;
            if uniqueAPLength
               templates{Ntemplates}   = sp;
            else
               templates(:,Ntemplates) = sp;
            end
            templateAPs(Ntemplates) = cell(1);
            templateAPs{Ntemplates} = sp;
            Nsamples(Ntemplates)    = 1;
            alignment_ind_templates(Ntemplates) = alignment_inds(i);

         % Get template the spike is maximally correlated with, add it to the
         % group & calculate the new mean. If matching on corr we don't care
         % about amplitude so normalise, else big APs will drown out small APs.
         % If matching on cov we do care about AP size so leave as is.
         else
            [~,k] = max(rho);
            if strcmpi( matchtype, 'corr' ) && normAPs
               sp = sp / std(sp); 
            end

            if uniqueAPLength
               [ templates{k}, templateAPs{k} ,alignment_ind_templates(k)] = addSpikeToTemplate( templates{k}, sp,alignment_ind_templates(k),alignment_inds(i), templateAPs{k} );
            else
               templates(:,k) = (templates(:,k)*Nsamples(k) + sp) / (Nsamples(k)+1);
               templateAPs{k}(:,(Nsamples(k)+1)) = sp;
            end
            Nsamples(k) = Nsamples(k) + 1;              
         end
         
          % We don't want to number of templates to get too long - it slows
          % down the code. every so often, 1000 spikes, we should cull
          % templates that have accrued almost no spikes
          TestPeriod = 100;
          if mod(i,TestPeriod)==0
            [templates,templateAPs,Nsamples,alignment_ind_templates,Ntemplates] = removeUncommonTemplates(templates,alignment_ind_templates,templateAPs,Nsamples,Ntemplates,10,0.01,TestPeriod);
          end
%           if uniqueAPLength
%             [Ntemplates,templates,templateAPs,Nsamples] = do_culling_for_speed(i,N,200,Ntemplates,templates,templateAPs,Nsamples);
%           end
      end
     
      
      
   end % end for each spike
   
   % get rid of templates with one or two spikes, as they're wacky & messing
   % up with our merging of templates because we end of correlating on a
   % very short piece of the spike
   tooSmall = Nsamples <= 2; 
   if uniqueAPLength
      templates( tooSmall )   = [];
      templateAPs( tooSmall ) = [];
      Nsamples( tooSmall )    = [];
      alignment_ind_templates( tooSmall) = [];
      Ntemplates              = Ntemplates - sum( tooSmall );

      % clean each template by removing time samples with very few spikes
      % contributing to it
      [templates, templateAPs,alignment_ind_templates] = cleanTemplates( templates, templateAPs,alignment_ind_templates, 2/3 );
      
   else
      templates(:, tooSmall ) = [];
      templateAPs( tooSmall ) = [];
      Nsamples( tooSmall )    = [];
      alignment_ind_templates(tooSmall) = [];
      Ntemplates              = Ntemplates - sum( tooSmall );
   end
   
   % now we've finished processing spikes, merge APtemplates that are
   % correlated more than the required threshold - AP templates may have
   % become more similar over time, as spikes have been added
   % obviously we change the correlation every time we merge APs, so
   % randomise the order we merge in, & only merge once
   cprintf( 'keywords', '\tMerging similar family templates...\n');
   [ templates, templateAPs, alignment_ind_templates,Nsamples ] = mergeSimilarTemplates( templates, templateAPs, Nsamples, matchthresh, uniqueAPLength,alignment_ind_templates );
   
   if uniqueAPLength
      % do a final clean of the templates, removing samples with too few
      % contributors to the ensemble
      [templates, templateAPs,alignment_ind_templates] = cleanTemplates( templates, templateAPs,alignment_ind_templates, 2/3 );
   end
   
   % now order AP families by how many APs went into them
   [~, order]   = sort(Nsamples, 'descend');
   if uniqueAPLength
      templates = templates(order);
   else
      templates = templates(:, order);
   end
   templateAPs = templateAPs(order);
   Nsamples     = Nsamples(order);
   alignment_ind_templates = alignment_ind_templates(order);
   cprintf( 'keywords', '\tFinished extracting templates\n');
end

%Based on some probabilistic process, we should remove APs that havent got
%enough spikes - statistically improbable
function [templates,templateAPs,Nsamples,alignment_ind_templates,Ntemplates] = removeUncommonTemplates(templates,alignment_ind_templates,templateAPs,Nsamples,Ntemplates,spike_thresh,min_prop_firing,TestPeriod)
    
    % Calculate the probability of a template being culled based on the
    % size of the spike threshold, and the minimum proportion of spikes we
    % expect from each unit. Templates will be tested against that
    % threshold with this probability. Thus, the longer a template is
    % around, the more likely it is to be tested, but also, the more likely
    % it is to have accumulated enough spikes to have passed the test.
    prob = TestPeriod*min_prop_firing/spike_thresh;
    to_remove = Nsamples<spike_thresh & rand(size(Nsamples))<prob;
    templates(to_remove) = [];
    templateAPs(to_remove) = [];
    Nsamples(to_remove) = [];
    alignment_ind_templates(to_remove) = [];
    Ntemplates = Ntemplates-sum(to_remove);
   

end

% Remove any part of a template where there's only 1 or a few spikes 
% contributing to that part of the template - assume this is noise
function [templates, templateAPs,alignment_ind_templates] = cleanTemplates( templates, templateAPs,alignment_ind_templates, perc )
   % can call this function with a single template
   if ~iscell( templates )
      singleTemp  = true;
      templates   = { templates };
      templateAPs = { templateAPs };
   else
      singleTemp  = false;
   end
   nTemp = length( templates );
   for ti=1:nTemp
      % get min number of spikes at each time sample based on min 
      % percentage of spikes required at each point
      nAPs  = size( templateAPs{ti}, 2 );
      Nreqd = floor( nAPs * perc );
      if nAPs > Nreqd
         Nensemble = sum( ~isnan( templateAPs{ti} ), 2 );
         templates{ti}( Nensemble < Nreqd )      = [];
         templateAPs{ti}( Nensemble < Nreqd, : ) = [];
         %need to shift alignment ind if we're cutting stuff off the start
         alignment_ind_templates(ti) = alignment_ind_templates(ti) - find((Nensemble < Nreqd)==0,1)+1;
      end
   end
   
   if singleTemp
      templates   = templates{:};
      templateAPs = templateAPs{:};
   end
end

% merge AP templates that are highly correlated
% component APs for each template are already centred around the peak
function [mergedTemp, mergedComp,alignment_inds_templates] = mergeTemplates( templates, componentAPs,alignment_inds_templates)
   nTemp    = length( templates );
   if nargin<3
      peakind  = getMaxInd( templates(:) ); % find index of template peak 
   else
      peakind =  alignment_inds_templates;
      %need to calculate a single alignment index for the new merged set,
      %since we're prepending nans, this will be the maximal alignment ind"
      alignment_inds_templates = max( peakind(:) );
   end
   nT       = cellfun(@length, templates(:) );    % length of template
   tbefore  = max( peakind(:) );
   tafter   = max( nT(:) - peakind(:) );
   
   if tbefore == 0 || tafter == 0
      str = sprintf( 'identifyUniqueTemplates:mergeTemplates:error: tbefore=%d, tafter=%d\n', tbefore, tafter );
      cprintf( 'Keywords*', str );
   end
   
   % append NaNs to template spikes or new spike to make them the same size
   prependfn= @( sp, N ) [ NaN( N, size(sp,2)); sp ];
   appendfn = @( sp, N ) [ sp; NaN( N, size(sp,2)) ];
   
   try
      tempmean = zeros( tbefore + tafter, nTemp );
   catch ME
      str = sprintf( 'identifyUniqueTemplates:mergeTemplates:error: tbefore=%d, tafter=%d\n', tbefore, tafter );
      cprintf( 'Keywords*', str );
      runtimeErrorHandler( ME );
   end
   
   try
      tmpAPs = componentAPs; 
      %Align all the templates
      for i=1:nTemp
         % if this template has fewer samples before peak, prepend NaNs 
         if peakind(i) < tbefore  
            componentAPs{i} = prependfn( componentAPs{i}, tbefore - peakind(i) );
         end
         if nT(i) - peakind(i) < tafter
         % if this template has fewer samples after peak, append NaNs 
            componentAPs{i} = appendfn( componentAPs{i}, abs( nT(i) - peakind(i) - tafter ) );
         end
         if size( componentAPs{i}, 1 ) ~= ( tafter + tbefore )
            str = sprintf( 'identifyUniqueTemplates:mergeTemplates:error: template is wrong size, ignoring...(%d instead of %d)\n', ...
                           size( componentAPs{i}, 1 ), ( tafter + tbefore ) );
            cprintf( 'Keywords*', str );
         end
         tempmean(:,i) = nanmean( componentAPs{i}, 2 );
        
      end
      
      try
          %merge the component APs then generate new template
         [mergedTemp,mergedComp] = CalculateTemplateFromAPs([ componentAPs{:} ],0.6,0.05,0,[]);
        
%          % mergedComp = horzcat( componentAPs{:} );
%          mergedComp = [ componentAPs{:} ]; 
% 
%          % now take an ensemble average of template, where the number of samples
%          % contributing to each time point may change across the spike, depending
%          % on the length of the constituent spikes 
%          mergedTemp = nanmean( mergedComp, 2 ); 
% 
%          perc  = 0.3;
%          nAPs  = size( mergedComp, 2 );
%          Nreqd = floor( nAPs * perc );
%          if nAPs > Nreqd
%             Nensemble = sum( ~isnan( mergedComp ), 2 );
%             mergedTemp( Nensemble < Nreqd )    = [];
%             mergedComp( Nensemble < Nreqd, : ) = [];
%          end
         
      catch M
         % mean that doesn't ditch samples with few inputs from the components
         N   = cellfun(@(t) size( t, 2 ), componentAPs );
         tmp = sum( tempmean .* repmat(N(:)', [tbefore+tafter 1]) / sum(N), 2);
         str = sprintf( 'ignoring number of samples in the ensemble when calculating template avg\n' );
         str = getCatchMEstring( M, str );
         cprintf( 'Keywords*', str );
      end
      
   catch ME
      str = sprintf( 'identifyUniqueTemplates:mergeTemplates:error calculating mean merged template\n' );
      str = getCatchMEstring( ME, str );
      cprintf( 'Keywords*', str );
      runtimeErrorHandler( ME );
   end
end


function [ templates, templateAPs,alignment_ind_templates, Nsamples ] = mergeSimilarTemplates( templates, templateAPs, Nsamples, matchthresh, uniqueAPLength,alignment_ind_templates)
   if uniqueAPLength
      Ntemp = length( templates );
   else
      Ntemp = size( templates, 2 );
      nT    = size( templates, 1 );
   end

   try
       
      %Before we do merging, we should align all the templates at their
      %peaks
      if uniqueAPLength
          %calculate full extent of template (not just parts with 2/3 samples
          full_templates = cellfun(@(i) {nanmean(i,2)},templateAPs);
          %weight indices of these templates based on proportion of samples
          template_weights =  cellfun(@(i) {sum(~isnan(i),2)/size(i,2)},templateAPs);
          %find peak location in template
          nT         = cellfun( @length, templates ); % num samples in each template
          if nargin<6
            peakind    = getMaxInd( templates, 1 );
          else
              peakind = alignment_ind_templates;
             
          end
          tafter     = ( nT - peakind );  
          max_after = max(tafter);
          max_before = max(peakind);
          %pad each template with nans so they are the same length
          aligned_templates = arrayfun(@(i) {[zeros(max_before-peakind(i),1);full_templates{i};zeros(max_after-tafter(i),1)]},1:length(full_templates));
          aligned_templates = [aligned_templates{:}];


          % make a matrix of weights for each index of each template (based on
          % the proportion of samples that have data at that index - this will
          % be used for weighting the correlation between the two templates,
          % and for combining the templates together
          aligned_weights = arrayfun(@(i) {[zeros(max_before-peakind(i),1);template_weights{i};zeros(max_after-tafter(i),1)]},1:length(full_templates));
          aligned_weights = [aligned_weights{:}];
      else
          %don't need to align templates because already aligned
          aligned_templates = templates;
          aligned_weights = ones(size(templates));
      end
      
      %merge
      merge_setting = 'pairwise_ordered';%'pairwise_ordered'
      switch merge_setting
          case 'random'
              %Random order merge 
              nAPs  = Ntemp; % start with all templates, and reduce as we merge
              order = randperm(nAPs, nAPs); % randomise merge order

              %method for merging:
              %do weighted pairwise correlation of aligned templates (where weights relate to
              % proportion of samples with data). for those that exceed threshold,
              % merge templates by compining samples and recalculating
              for i=3:nAPs
                 if i>nAPs, break; end % reducing number of APs each time through so check
                 jj      = order(i);

                 %Step 1: find templates to merge based on correlation of templates
                 if uniqueAPLength
                     rho = corr_weighted(aligned_weights,aligned_templates);
                     merge = find(rho(jj,:)>matchthresh);
                 else
                     rho     = corr(  templates(:,jj), templates );
                     merge   = find( rho > matchthresh );
                 end


                 %If there are templates to merge then we need to combine them and
                 %remove their old indices
                 if numel(merge) > 1 

                     %Step 2: recalculated aligned templates for the merged template as
                     %well as recalculating templates and templateAPs
                     if uniqueAPLength
                        %recalculate templates as weighted sum of original templates where
                        %weights are number of samples with data in that spot
                        aligned_templates(:,jj) = sum(aligned_templates(:,merge).*aligned_weights(:,merge).*Nsamples(merge),2)./sum(Nsamples(merge));
                        aligned_weights(:,jj) = sum(aligned_weights(:,merge).*Nsamples(merge),2)./sum(Nsamples(merge));
                        [templates{jj}, templateAPs{jj}] = mergeTemplates( templates( merge ), templateAPs( merge ) );
                     else
                         N     = Nsamples(merge);
                         aligned_templates(:, jj) = sum( aligned_templates(:, merge) .* repmat(N(:)', [nT 1]) / sum(N), 2); 
                         templates(:, jj) = aligned_templates(:, jj);
                         templateAPs{jj}  = cell2mat( templateAPs(merge) ); 
                     end

                     %step 3: remove elements that were merged together and adjust
                     %order to reflect that
                     ditch = merge;
                     ditch(ditch==jj) = [];
                     aligned_templates(:,ditch) = [];
                     aligned_weights(:,ditch) = [];
                     Nsamples(jj)        = sum( Nsamples(merge) ); % add samples from all matching templates      
                     Nsamples(ditch)     = [];    % get rid of merged template num samples
                     nAPs                = nAPs - numel(ditch);
                     if uniqueAPLength
                        templates(ditch) = [];
                        templateAPs(ditch) = [];
                     else
                        templates(:, ditch) = []; % get rid of merge AP templates
                        templateAPs(ditch)  = [];
                     end
                     order               = arrayfun(@(ii) ternaryOp(any(ii==ditch), 0, ii), order);
                     order(order==0)     = [];
                      % udpate order indices by counting how many of the merged
                      % templates were sitting before each template
                     order               = arrayfun(@(ii) ii - sum(ii > ditch), order);
                 end 
              end
          case 'pairwise_ordered'
                %Random order merge 
              nAPs  = Ntemp; % start with all templates, and reduce as we merge
              order = randperm(nAPs, nAPs); % randomise merge order

              %method for merging:
              %do weighted pairwise correlation of aligned templates (where weights relate to
              % proportion of samples with data). for those that exceed threshold,
              % merge templates by compining samples and recalculating
           
              if uniqueAPLength
                   rho = corr_weighted(aligned_weights,aligned_templates);
                   max_rho = max(rho.*~eye(size(rho)),[],'all');
              else
                   rho     = corr(  templates, templates );
                   max_rho = max(rho.*~eye(size(rho)),[],'all');
              end
              
              
              while max_rho>matchthresh
                 
                [merge(1),merge(2)] = find(rho==max_rho,1);
                 jj = merge(1);
                     %Step 2: recalculated aligned templates for the merged template as
                     %well as recalculating templates and templateAPs
                     if uniqueAPLength
                        %recalculate templates as weighted sum of original templates where
                        %weights are number of samples with data in that spot
                        aligned_templates(:,jj) = sum(aligned_templates(:,merge).*aligned_weights(:,merge).*Nsamples(merge),2)./sum(Nsamples(merge));
                        aligned_weights(:,jj) = sum(aligned_weights(:,merge).*Nsamples(merge),2)./sum(Nsamples(merge));
                        [templates{jj}, templateAPs{jj},alignment_ind_templates(jj)] = mergeTemplates( templates( merge ), templateAPs( merge ),alignment_ind_templates(merge) );
                     else
                         N     = Nsamples(merge);
                         aligned_templates(:, jj) = sum( aligned_templates(:, merge) .* repmat(N(:)', [nT 1]) / sum(N), 2); 
                         templates(:, jj) = aligned_templates(:, jj);
                         templateAPs{jj}  = cell2mat( templateAPs(merge) ); 
                     end

                     %step 3: remove elements that were merged together and adjust
                     %order to reflect that
                     ditch = merge;
                     ditch(ditch==jj) = [];
                     aligned_templates(:,ditch) = [];
                     aligned_weights(:,ditch) = [];
                     Nsamples(jj)        = sum( Nsamples(merge) ); % add samples from all matching templates      
                     Nsamples(ditch)     = [];    % get rid of merged template num samples
                     nAPs                = nAPs - numel(ditch);
                     if uniqueAPLength
                        templates(ditch) = [];
                        templateAPs(ditch) = [];
                        alignment_ind_templates(ditch) = [];
                     else
                        templates(:, ditch) = []; % get rid of merge AP templates
                        templateAPs(ditch)  = [];
                     end
                    
                  
                   if uniqueAPLength
                       rho = corr_weighted(aligned_weights,aligned_templates);
                       max_rho = max(rho.*~eye(size(rho)),[],'all');
                   else
                       rho     = corr(  templates, templates );
                       max_rho = max(rho.*~eye(size(rho)),[],'all');
                   end
                 
              
              
              end 
              
   
      end
   catch ME
      str = getCatchMEstring( ME, 'end: ' );
      runtimeErrorHandler(ME);
   end
end




% Extract correlation only of templates similar to current template, both
% in number of samples before & after the peak, and in shape. 
% (necessary because can only correlate on the same length spikes, so if we
% extract out common samples before and after the peak, we sometimes end up
% correlating on just a short straight line, so everything is collapsed
% into a single template)
function merge = getSimilarTemplates( templates, jj, uniqueAPLength, matchthresh,alignment_ind_templates )
   range = 3; % allow this number of samples different btwn samples before & after peak 
   try
      % if templates are all the same length simply correlate & threshold
      if ~uniqueAPLength
         rho     = corr(  templates(:,jj), templates );
         merge   = find( rho > matchthresh );
         return;
      end
      
      % if templates are all different lengths identify the templates with
      % similar number of samples before & after the peak, & then correlate
      nT         = cellfun( @length, templates ); % num samples in each template
      if nargin<5
        peakind  = getMaxInd( templates, 1 );
      else
         peakind = alignment_ind_templates;
      end
      tafter     = ( nT - peakind );      % number of samples after peak
      curr_temp  = templates{jj};         % current template
      curr_nT    = length( curr_temp );   % length of current template
      curr_peak  = peakind(jj);           % peak index of current template
      curr_tafter= tafter(jj);            % num samples after peak for curr template 
      similar    = find( ( (curr_peak  -range) < peakind & peakind < (curr_peak+range) ...
                        &  (curr_tafter-range) <  tafter &  tafter < (curr_tafter+range) ) ); 
      tbefore    = min( peakind(similar) );
      tafter     = min( nT(similar) - peakind(similar) );
      temps      = arrayfun( @(i) templates{i}( peakind(i)-tbefore+1 : peakind(i)+tafter ), similar, 'uni', false ); 
      temps      = cell2mat( temps );

      % gotta calc corr or cov separately for each template because
      % they're different lengths
      temps     = arrayfun( @(i) templates{i}( peakind(i)-tbefore+1 : peakind(i)+tafter ), ...
                            similar, 'uni', false );
      curr_temp = curr_temp( curr_peak-tbefore+1 : curr_peak+tafter );
      % make matrix copy of templates using min duration before & after peak
      temps     = cell2mat( temps ); 
      rho       = corr( curr_temp, temps );
      merge     = rho > matchthresh; % threshold correlation
      merge     = similar( merge );  % gotta convert back to orig template indices

   catch ME
      str = getCatchMEstring( ME, 'I hate this runtime error business\n' );
   end
end














