function path = getFilePath(varargin)
   if nargin >= 1
      option = varargin{1};
   else
      option = 'not';
   end
   
   a = which('SpikeExtractionTool');
   locs = strfind(a, filesep);
   path = a(1:locs(end));
   [~,att] = fileattrib(path);

   if strcmp('log', option)
      if ~exist(path,'file') || ~att.UserWrite
         path = userpath;
         [~,att] = fileattrib(path);
         if ~exist(path,'file') || ~att.UserWrite
            path = ctfroot;
            [~,att] = fileattrib(path);
            if ~exist(path,'file') || ~att.UserWrite
               displayErrorMsg('Couldn''t find a folder to save a Log File. Please chose one.');
               path = uigetdir();
            end
         end
      end
   end
end