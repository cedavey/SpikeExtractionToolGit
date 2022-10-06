% COMPARE_TEMPLATES Find the closest match of templates with current spike.
% Not the best correlation/covariance, but given 2 similarly matching
% rho's, chose the one with a smaller size difference.
%
% Artemio Soto-Breceda | 22-August-2019
function [k, d] = compare_templates(rho, curr_spike, APtemplates, APfamilies, matchthresh)
   d = [];
   if iscell( APtemplates ) 
      uniqueAPLength = true;
   else
      uniqueAPLength = false;
   end


   % get difference btwn peaks and difference btwn troughs
   compare  = @(s1, s2) 2*( abs( max(s1) - max(s2) ) ) + abs( min(s1) - min(s2) ); % wgt peak more impt
   compare  = @(s1, s2) abs( max(s1) - max(s2) ) + abs( min(s1) - min(s2) ); % wgt peak & trough the same
   rho_     = rho/max(rho); % normalize rho to current largest

   nfam     = cellfun( @length, APfamilies );

   % Check the templates with active families & within matchthresh % similarity
   similars = find(rho_>matchthresh & nfam>0);
   if numel( similars ) <= 1
      if isempty(similars)
         [~,k] = max(rho);
      else
         k = similars;
      end
      return;
   end

   % only consider templates that are similar to max rho
   if uniqueAPLength
      APtemplates = APtemplates(similars);
      nap = length( APtemplates );
   else
      APtemplates = APtemplates(:,similars);
      nap = size( APtemplates, 2 );
   end
   APfamilies = APfamilies( similars );
   rho_red    = rho(similars);   % rho just for templates with similar corr
   [~, k_red] = max(rho_red);    % template index incl. only templates with similar corr

   d    = inf( nap, max(nfam) ); % distance between peaks & troughs for each axon
   
   for ai = 1:nap

      % get smallest correlation with families that match this template
      for fi=1:length( APfamilies{ai} )
         fspike = APfamilies{ai}{fi}.meanspike;
         
         % curr_spike & fspike may hae diff lengths
         if uniqueAPLength 
            % if spikes & templates can have diff lengths, use the smaller length
            % to calculate the match, but centre both around the peak 
            [ curr_spike, fspike ] = alignDiffLengthSpikes( curr_spike, fspike );

            % gotta calc corr or cov separately for each template because
            % they're different lengths
            try
               d(ai,fi) = compare( curr_spike, fspike );
            catch
               % something has gone wrong - return largest corr template
               [~, k] = max(rho);
               return
            end

         else
            d(ai,fi) = compare( curr_spike, fspike );
         end
         
      end
   end

   % if optimal rho has distance almost as good as others, don't change k
   d_temps    = min(d,[],2);     % get best distance for each template
   d_currtemp = d_temps(k_red);  % best distance for template with best corr
   d_others   = d_temps;         % to contain distance btwn non-optimal templates & spike
   d_others(k_red) = [];         % get rid of current optimal template
   % if distance to template with best corr is only marginally worse than
   % the others, prioritise the better correlation rather than the better
   % distance between peaks and troughs
   rho_withfam= rho(similars);   % consider corr of templates with family
   [~, k_rho] = max(rho_withfam);% original optimal template, based only on corr
   k_rho      = similars(k_rho); % get k out of templates with good corr & families
   
   % template with min distance differs from optimal rho template, but the
   % difference is marginal --> revert to prioritising correlation 
   if min(d_temps)~=d_currtemp && all(d_currtemp < d_others*1.05)
      k = k_rho;
      return;
   end

   % calculate best corr, but this is for template & family within the
   % template, so you need to make k just a function of template, not family
   [~, k_rho] = max(rho);        % original optimal template, based only on corr
   [~, k]     = min( d(:) );     % best distance across templates & families
   [k, fi]    = ind2sub( [nap max(nfam)], k); % get best template, ignore family
   k          = similars(k);     % index including all templates, not just similar templates
   if k~=k_rho
      str = sprintf( 'compare_templates: changing k from %d to %d\n', k_rho, k );
%       cprintf( 'Keywords*', str );
   end
end

% Separate a spike into positive and negative parts
function [vpos, vneg] = divide_templates(v)   
   for ii = 1:size(v,2)
      [~, k] = max(v(:,ii));
      d = diff(v(:,ii) > 0);
      try
         i = find(d(1:k)   == 1,1, 'last');
         j = find(d(k:end) == 1,1, 'first') + k;
      catch E
         if strcmp('MATLAB:badsubscript', E.identifier)
            k = numel(d) - 1;
            i = find(d(1:k)   == 1,1, 'last');
            j = find(d(k:end) == 1,1, 'first') + k;
         end
      end
      v(1:i, ii)   = 0;
      v(j:end, ii) = 0;
   end
   
   vpos = v;
   vneg = v;
   
   vpos(v < 0) = 0;
   vneg(v > 0) = 0;
   vneg = abs(vneg);
end



