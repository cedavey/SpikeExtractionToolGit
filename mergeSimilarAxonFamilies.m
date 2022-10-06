% Check a new axon for similarity with existing axons, to ensure that it
% wasn't created from a single stray axon. If it is similar we need to
% merge it with the family it's similar to, and then update variables that
% track family info
function [APfamilies, is_new_axon, num_new_axons] = mergeSimilarAxonFamilies( APfamilies, this_axon, is_new_axon, num_new_axons, ...
                                                    rhothresh, kappa_pos, kappa_neg )
   k    = this_axon(1); 
   a    = this_axon(2);
   axon = APfamilies{k}{a}; 

   Nfam = length(APfamilies);

   % for through each axon family & check its similarity with this axon
   % we shouldn't ever merge & then find that the merged family is similar
   % to another family, as they would have been merged already (the number
   % of spikes in the axon family being checked is small, so shouldn't
   % change the mean too much)
try
   for ki=1:Nfam
      Naxon = length(APfamilies{ki});
      for ai=1:Naxon
         % don't want to compare with itself! or with other new axons
         if ~(ki==k && ai==a) && ~is_new_axon(ki,ai) 
            fam = APfamilies{ki}{ai}; 
            [similar, mergedfam] = compareAxonFamilies( axon, fam, rhothresh, kappa_pos, kappa_neg );
            if similar
               APfamilies{ki}{ai} = mergedfam; 
               % if this is the only axon in the family, empty the family,
               % else shift the other axons down
               if length(APfamilies{k}) > a
                    is_new_axon(k,1:end-1) = [   is_new_axon(k,1:a-1)   is_new_axon(k,a+1:end) ];
                  num_new_axons(k,1:end-1) = [ num_new_axons(k,1:a-1) num_new_axons(k,a+1:end) ];
                    is_new_axon(k,end) = 0;
                  num_new_axons(k,end) = 0;
               else
                  % don't remove this axon field in case other families
                  % need the matrix this large, but set to 0
                    is_new_axon(k,a)   = 0;
                  num_new_axons(k,a)   = 0;
               end
               APfamilies{k}(a) = [];
            end
            return;
         end
      end

   end
catch ME
   disp('Problem with the merge');
end
   
end