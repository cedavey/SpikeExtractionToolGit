% Removes a string from the list of stat names (e.g. for if a mask or time
% series is removed). Assumes that where a string is repeated an index has
% been appended to the end, e.g. voltage4
% Inputs:
%   names   - cell array of currently used base names & frequency of use
%   rmstr   - string to remove from list
function names = removeStringFromList(names, rmstr)
   for i=1:size(names,1)
      str = names{i,1};
      if ~isempty(strfind(rmstr, str))
         % if only 1 use of str then remove from names
         if names{i,2}==1
             names(i,:) = [];
             
         % if more than 1 use of str then if rmstr is the last occurence
         % of it (i.e. has index = frequency of str) then reduce
         % frequency by 1
         elseif names{i,2}==str2num(rmstr(end)) % rmstr index is last occurence
            % if there's more than 1 use of str, & rmstr is an earlier
            % occurence of it, then just leave frequency as is, else if we
            % reduce it by 1 when the str next occurs we'll repeat the last index
             names{i,2} = names{i,2} - 1;
         end
         break;
      end
   end
end

