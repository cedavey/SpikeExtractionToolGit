% [APtemplates, Nsamples, componentAPs] = identifyUniqueTemplates(spikes, matchtype, matchthresh)
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
%  matchtype    - 'corr' or 'cov': use corr to match only on shape, cov to 
%                 match on shape & size (format: string, defaults" 'corr')
%  matchthresh  - correlation of spike with template must be less than this 
%                 to be considered unique & a new template initiated
%                (format: float scalar, default: 0.95)
%  normAPs      - if true normalise AP templates for future spike matching
%                (format: boolean, default: true)
%
% Outputs:
%  uniqueAPs    - copy of all unique action potential waveforms 
%                 (format: time x k, where k is the number of unique APs      )
%              OR (depending on whether the input 'spikes' is a cell or matrix)
%                 (format: cell x k, and time x 1 within each cell, but the   )
%                  duration of time changes for each template                 )
%  Nsamples     - number of APs used in construction of mean AP waveform (kx1)
%                 (cell array of matrices, 1 cell per template, and a
%                 single vector of the template, but templates have diff
%                 length vectors)
%  componentAPs - keeps a list of all spikes contributing to the templates 
%                 AP (cell array of matrices, 1 cell per template, and 
%                 matrix size template_length x num_spikes)
function [templates, Nsamples, templateAPs] = identifyUniqueTemplates(spikes, varargin)
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
   Nsamples       = 1;  % number of AP samples for each waveform
   Ntemplates     = 1;  % record the number of templates we have

   % get match function
   switch matchtype
      case 'corr'
         matchfn = @( sp, temp ) corr( temp, sp );

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
   for i=2:N
      sp  = getspike(i); 
      rho = spikeCorrWithTemplates( sp, templates, matchfn );
      
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
               [ templates{k}, templateAPs{k} ] = addSpikeToTemplate( templates{k}, sp, templateAPs{k} );
            else
               templates(:,k) = (templates(:,k)*Nsamples(k) + sp) / (Nsamples(k)+1);
               templateAPs{k}(:,(Nsamples(k)+1)) = sp;
            end
            Nsamples(k) = Nsamples(k) + 1;              
         end
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
      Ntemplates              = Ntemplates - sum( tooSmall );

      % clean each template by removing time samples with very few spikes
      % contributing to it
      [templates, templateAPs] = cleanTemplates( templates, templateAPs, 2/3 );
      
   else
      templates(:, tooSmall ) = [];
      templateAPs( tooSmall ) = [];
      Nsamples( tooSmall )    = [];
      Ntemplates              = Ntemplates - sum( tooSmall );
   end
   
   % now we've finished processing spikes, merge APtemplates that are
   % correlated more than the required threshold - AP templates may have
   % become more similar over time, as spikes have been added
   % obviously we change the correlation every time we merge APs, so
   % randomise the order we merge in, & only merge once
   cprintf( 'keywords', '\tMerging similar family templates...\n');
   [ templates, templateAPs, Nsamples ] = mergeSimilarTemplates( templates, templateAPs, Nsamples, matchthresh, uniqueAPLength );
   
   if uniqueAPLength
      % do a final clean of the templates, removing samples with too few
      % contributors to the ensemble
      [templates, templateAPs] = cleanTemplates( templates, templateAPs, 2/3 );
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
   cprintf( 'keywords', '\tFinished extracting templates\n');
end


% Remove any part of a template where there's only 1 or a few spikes 
% contributing to that part of the template - assume this is noise
function [templates, templateAPs] = cleanTemplates( templates, templateAPs, perc )
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
      end
   end
   
   if singleTemp
      templates   = templates{:};
      templateAPs = templateAPs{:};
   end
end

% merge AP templates that are highly correlated
% component APs for each template are already centred around the peak
function [mergedTemp, mergedComp] = mergeTemplates( templates, componentAPs )
   nTemp    = length( templates );
   peakind  = getMaxInd( templates(:) ); % find index of template peak 

   nT       = cellfun(@length, templates(:) );    % length of template
   tbefore  = max( peakind(:) );
   tafter   = max( nT(:) - peakind(:) );
   
   if tbefore == 0 || tafter == 0
      str = sprintf( 'identifyUniqueTemplates:mergeTemplates:error: tbefore=%d, tafter=%d\n', tbefore, tafter );
      cprintf( 'Keywords*', str );
      runtimeErrorHandler( ME );
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
         % mergedComp = horzcat( componentAPs{:} );
         mergedComp = [ componentAPs{:} ]; 

         % now take an ensemble average of template, where the number of samples
         % contributing to each time point may change across the spike, depending
         % on the length of the constituent spikes 
         mergedTemp = nanmean( mergedComp, 2 ); 

         perc  = 0.3;
         nAPs  = size( mergedComp, 2 );
         Nreqd = floor( nAPs * perc );
         if nAPs > Nreqd
            Nensemble = sum( ~isnan( mergedComp ), 2 );
            mergedTemp( Nensemble < Nreqd )    = [];
            mergedComp( Nensemble < Nreqd, : ) = [];
         end
         
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


function [ templates, templateAPs, Nsamples ] = mergeSimilarTemplates( templates, templateAPs, Nsamples, matchthresh, uniqueAPLength )
   if uniqueAPLength
      Ntemp = length( templates );
   else
      Ntemp = size( templates, 2 );
      temps = templates; % gotta make a copy cuz if APtemplates is a cell we made a matrix version
      nT    = size( templates, 1 );
   end

   try
      nAPs  = Ntemp; % start with all templates, and reduce as we merge
      order = randperm(nAPs, nAPs); % randomise merge order
      for i=3:nAPs
         if i>nAPs, break; end % reducing number of APs each time through so check
         jj      = order(i);
         merge   = getSimilarTemplates( templates, jj, uniqueAPLength, matchthresh );

         if numel(merge) > 1 % number of similar templates must be more than current template
            ditch = merge; % don't wanna ditch the template we merged into
            ditch( merge==jj ) = []; % don't merge, & therefore ditch, current temp (merge INto it!)
            % when calculating merged AP template, weight avg by num samples
            N     = Nsamples(merge); % weight for each template
            % if AP templates are all same length & in a matrix, update & 
            % merge differently to if they're different lengths & kept in cells
            if uniqueAPLength
               % if templates have unique constituent AP length, merge from
               % component spikes & take ensemble average
               [templates{jj}, templateAPs{jj}] = mergeTemplates( templates( merge ), templateAPs( merge ) );
               templates(ditch) = []; % get rid of merge AP templates
            else
               temps(:, jj)     = sum( temps(:, merge) .* repmat(N(:)', [nT 1]) / sum(N), 2); 
               templates(:, jj) = temps(:, jj);
               templates(:, ditch) = []; % get rid of merge AP templates
               templateAPs{jj}  = cell2mat( templateAPs(merge) ); 
               temps(:, ditch)  = [];    % get rid of merged templates
            end
            templateAPs(ditch)  = [];
            Nsamples(jj)        = sum( Nsamples(merge) ); % add samples from all matching templates      
            Nsamples(ditch)     = [];    % get rid of merged template num samples
            nAPs                = nAPs - numel(ditch);
            % remove indices of merged AP templates 
            mergeind            = find(ditch);
            order               = arrayfun(@(ii) ternaryOp(any(ii==ditch), 0, ii), order);
            order(order==0)     = [];
             % udpate order indices by counting how many of the merged
             % templates were sitting before each template
            order               = arrayfun(@(ii) ii - sum(ii > ditch), order);
         end
      end

   catch ME
      str = getCatchMEstring( ME, 'fucking runtime error handler !!\n' );
   end
end

% Extract correlation only of templates similar to current template, both
% in number of samples before & after the peak, and in shape. 
% (necessary because can only correlate on the same length spikes, so if we
% extract out common samples before and after the peak, we sometimes end up
% correlating on just a short straight line, so everything is collapsed
% into a single template)
function merge = getSimilarTemplates( templates, jj, uniqueAPLength, matchthresh )
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
      peakind    = getMaxInd( templates, 1 );
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
      str = getCatchMEstring( ME, 'I fucking hate this runtime error business\n' );
   end
end

















