% COMPARE_TEMPLATES Find the closest match of templates with current spike.
% Not the best correlation/covariance, but given 2 similarly matching
% rho's, chose the one with a smaller size difference.
%
% Artemio Soto-Breceda | 22-August-2019
function k = compare_templates(rho, curr_spike, APtemplates, APfamilies)
   compare = @(s1, s2) 2*( abs( max(s1) - max(s2) ) ) + abs( min(s1) - min(s2) );
   rho_ = rho/max(rho); % normalize rho

   % Check the templates with 5% similarity
   similars = find(rho_ > 0.95);
   if numel( similars ) <= 1
      [~, k] = max(rho_);
      return;
   end

   APtemplates = APtemplates(:,similars);
   APfamilies  = APfamilies( similars );
   
   nap  = size( APtemplates, 2 );
   nfam = cellfun( @length, APfamilies );
   nosp = find( nfam == 0 ); % record templates with no families 
   nfam( nosp ) = 1; % make templates with no families have 1 family/spike, given by template
   d    = inf( nap, max(nfam) );
   
   % spk         = curr_spike/max(abs(curr_spike)); % normalize size  
   % [spk_pos, spk_neg] = divide_templates( spk );

   for ai = 1:size(APtemplates,2)
      if any( ai==nosp )
         APfamilies{ai}{1}.meanspike = APtemplates(:,ai);
      end
      for fi=1:length( APfamilies{ai} )
         fspike = APfamilies{ai}{fi}.meanspike;
         d(ai,fi) = compare( curr_spike, fspike );
         
%          APtemplates(:,ai) = APtemplates(:,ai) - mean( APtemplates(:,ai) );
%          APtemplates(:,ai) = APtemplates(:,ai) / max( abs(APtemplates(:,ai)) ); % Normalize size
      end
   end      
%    [tmps_pos, tmps_neg] = divide_templates( APtemplates );

   % weight peak as more important than trough
%    difference = 2 * ( max(tmps_pos) - max(spk_pos) ) + max(tmps_neg) - max(spk_neg);
%    [~, k] = min( abs(difference) );
%    k = similars(k);

   [~, k] = min( d(:) );
   [k,~] = ind2sub( [nap max(nfam)], k );
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



