% COMPARE_TEMPLATES Find the closest match of templates with current spike.
% Not the best correlation/covariance, but given 2 similarly matching
% rho's, chose the one with a smaller size difference.
%
% Artemio Soto-Breceda | 22-August-2019
function k = compare_templates(rho, curr_spike, tmps)

   rho_ = rho/max(rho); % Normalize rho.

   % Check the templates with 10% similarity
   similars = find(rho_ > 0.95);
   if numel(similars) > 1
      tmps = tmps(:,similars);
      spk = curr_spike/max(abs(curr_spike)); % Normalize size
      [spk_pos, spk_neg] = divide_templates(spk);
      
      for i = 1:size(tmps,2)
         tmps(:,i) = tmps(:,i) - mean(tmps(:,i));
         tmps(:,i) = tmps(:,i)/max(abs(tmps(:,i))); % Normalize size
      end      
      [tmps_pos, tmps_neg] = divide_templates(tmps);
      
      % difference = (sum(tmps_pos) - sum(spk_pos)) + (sum(tmps_neg) - sum(spk_neg));
      difference = max(tmps_pos) - max(spk_pos) + max(tmps_neg) - max(spk_neg);
      [~, k] = min(abs(difference));
      k = similars(k);
   else
      [~, k] = max(rho_);
   end

end

% Separate a spike into positive and negative parts
function [vpos, vneg] = divide_templates(v)   
   for ii = 1:size(v,2)
      [~, k] = max(v(:,ii));
      d = diff(v(:,ii) > 0);
      i = find(d(1:k) == 1,1, 'last');
      j = find(d(k:end) == 1,1, 'first') + k;
      v(1:i, ii) = 0;
      v(j:end, ii) = 0;
   end
   
   vpos = v;
   vneg = v;
   
   vpos(v < 0) = 0;
   vneg(v > 0) = 0;
   vneg = abs(vneg);
end