% [new_str, str_list] = check4RepeatString(new_str, str_list)
%
% (Function made initially for gui so it's written in this context)
% Makes sure there are no repeat strings in the cell array 'stats_str',
% because if there are the popup menu can't distinguish between those stats
% with the same name (weird huh?!). If there is a repeat then add an index
% to it. If you've added an index to a basename
% already, then you need to make sure you increase the index the next time,
% else the 2 indexed strings will be the same.
% Inputs:
%   str    - new unique stat name(s) - (basename_# where nec) 
%   names  - list of used basenames, & frequency of use
% Outputs:
%   str    - list of current unique names
%   names  - list of currently used basenames & frequency of use

% FIXME: rather than comparing all strings against each other every time,
% you should include the option of giving an index indicating where the new
% strings start, and then check only from there on. Save time.
function [new_str,str_list] = check4RepeatString(new_str, str_list)
    % format of the names cell array is:
    %   col1 = basename (string)
    %   col2 = current index
    % E.g. names{1,1} = 'test', names{1,2} = 1 --> there's already 1 stat
    % called 'test'
    
    % function written to allow addition of multiple string simultaneously
    % --> expects cell array. For single string, convert to cell
    if ischar(new_str), new_str = {new_str}; end
    nnames = size(str_list,1);
    % get number of strings we know about already
    if isempty(str_list)
       n = 0;
    else
       n = sum(cell2mat(str_list(:,2)));
    end
    nstats = length(new_str);
%     for i=n+1:nstats % start from new strings
    for i=1:nstats % start from new strings
        new_str = new_str{i};
        
        % if names{j} is found in str then we need to check where it ends
        % in str cuz, e.g., maybe there's arma & ar, or something, so ar is
        % found in arma, but it's not the same
        j=1;
        while j<=nnames
            base = str_list{j,1};
            k  = strfind(new_str,base); % start index of match
            k_ = k+length(base)-1;  % end index of match
            % if match is real then there's only an index number or nothing
            % (if this is 1st double-up) at the end
            if ~isempty(k) && k(end)==1 && length(new_str)==k_(1) % then they're exactly the same so match=yes
                str_list{j,2} = str_list{j,2}+1;
                new_str = [new_str '_v' num2str(str_list{j,2})];
                break;
            end
            j = j+1;
        end
        if j>nnames % no match so add str to basename list
            str_list{j,1}=new_str;
            str_list{j,2}=1;
            nnames = nnames+1;
        end
    end
return

% This version is trying to only use stats_str (not names), but I realised
% that I'd need to find out where the new set of strings started because
% I'm starting the search from the end, but if there's 2 new strings that
% match then their versions will be back-to-front (earlier string in list
% will have higher version number). This is a problem because I'm using the
% first occurence of a version number (going backwards) to know how many
% versions I have. But if an earlier string might have a higher version,
% then I'd end up giving the new string the wrong version.
%         while j>=1
%             str2 = stats_str{j};
%             % gotta search for 'str' pattern in str2, because there may
%             % exist either versioned copies of str already, or an exact
%             % match (if no duplicates until now), so str2 will be either
%             % shorter or equal in length to it's duplicates
%             % if we've already got multiple occurences of a string, they'll
%             % be versioned and therefore won't match.
%             k  = strfind(str2,str); % start index of match
%             k_ = k+length(str2)-1;  % end index of match
%             % if match is real then end index will be at the end of the
%             % strings (1st duplicate) or just before versioning info
%             % (duplicates exist already )
%             % e.g. str='ar', we may already have 'ar_v1' & 'ar_v2', 
%             % or we may only have 'ar'
%             if length(str)==k_ % then they're exactly the same so match
%                 match = 1;
%                 break;
%             elseif (k_==length(str)-3) && (all(str(k_+1:k_+2)=='_v') )
%                 match = 1;
%                 break;
%             end
%             j = j-1;
%         end
