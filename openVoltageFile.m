% Opens the 'open file' dialogue to get the file name from the user, then
% attempts to open the selected file & retrieve the requested size
% matrices (either 3D or 4D since we're working with image files)
% Inputs:
%   data_struct     - current user data structure
%   sz              - size of matrices considered valid
% Outputs:
%   data_struct     - updated data struct (containing new image volumes)
%   success         - 1 if we successfully obtained new images
function [data, success] = openVoltageFile(data, varargin)
   % Check if optional input
   if nargin > 1
      isBatch = varargin{1};
      if strcmp(isBatch, 'batch')
         isBatch = true;
         path = varargin{2};
         file = varargin{3};
      else
         isBatch = false;         
      end
   else 
      isBatch = false;
   end
   
   success = 0;
   
   if ~isBatch
      ftypes = {'*.mat;*.smr;*.txt', 'Voltage files (*.mat, *.smr, *.txt)'; ...
        '*.mat', 'M-files (*.mat)'; ...
        '*.smr','Analyze image files (*.img)'; ...
        '*.txt','Text files (*.txt)'; ...
        '*.*',  'All Files (*.*)'};

      user_dir = ternaryOp( isdir( data.last_dir ), data.last_dir, userpath ); 
      [flist, pathname, filtindex] = uigetfile( ftypes, ...
        'Load voltage timeseries', user_dir, 'Multiselect', 'on');
      % if user cancelled action pathname should be 0 so just return
      if isequal(pathname,0) && isequal(flist,0)
        success = -1; % success==-1 --> user cancelled
        return
      end
      % remember dir of chosen file so can default to it next time
      data.last_dir = pathname; 
      
      % convert to cell cuz if user selects multiple files then flist is a
      % cell array - commonise the type so can deal with all cases
      if ~iscell(flist), flist={flist}; end

      % get current time series details
      num_tseries = data.num_tseries;
      tseries_str = data.tseries_str;
      tseries     = data.tseries;
      used_names  = data.used_names;
   else
      flist = file;
      pathname = path;
      
      data.last_dir = pathname; 
      % get current time series details
      num_tseries = data.num_tseries;
      tseries_str = data.tseries_str;
      tseries     = data.tseries;
      used_names  = data.used_names;
   end
   
   for fname=flist
      % convert fname as cell to string
      fname    = fname{1};
      fullname = fullfile(pathname, fname); % builds full file name from parts   

      try
         [newdata, time] = readSpikeData(fullname, [], isBatch);
      catch ME
         newdata = []; time = [];
      end
      if isempty(newdata)
         str = sprintf('SpikeExtractionTool - LoadError: %s does not contain valid data, ignoring...', fname);
         printMessage('off', '*Errors', str);
        return;
      end
      
      % if voltage is a struct assume it's a previously saved SET structure
      if isstruct(newdata)
         SETdata = newdata;
         newdata = SETdata.data;
         time    = SETdata.time;
         
      % data is from an smr file or similar - assume it's a voltage timeseries
      else
         % use filename as data name to avoid confusion with multiple tseries
         [~, file, ext]  = fileparts( fname ); % remove extension
         SETdata.type    = 'voltage';
         SETdata.name    = file; % 'voltage';
         SETdata.dt      = time(2) - time(1);
         SETdata.params  = [];
      end
         
      if strcmpi(SETdata.type, 'voltage') && length(time)>1e4
         if ~isBatch
            response = userConfirmation('Truncate voltage timeseries?','Voltage length');
            switch lower(response)
              case 'yes'
                 [newdata, time] = truncateVoltage(newdata, time);
              case 'no'
            end
         end
         SETdata.data = newdata; 
         SETdata.time = time;
      end

      %% TO DO: allow saving & loading of voltage from gui with extra info
      try
         [name, used_names]  = check4RepeatString(SETdata.name, used_names);
         num_tseries         = num_tseries + 1;
         tseries_str{end+1}  = name;
         tseries{end+1}      = SETdata;
         [tlim, vlim]        = getTimeAndVoltageLimits(SETdata);
      catch ME
         return;
      end
   end

   data.num_tseries  = num_tseries;
   data.tlim         = tlim;
   data.vlim         = vlim;
   data.tseries_str  = tseries_str;
   data.tseries      = tseries;
   data.used_names   = used_names;
   data.last_tseries = data.curr_tseries; % shift curr tseries index to last tseries
   data.curr_tseries = num_tseries; % update curr time series index to new tseries
   data.last_dir     = pathname;    % update last directory selected by user
   
   % Check that none tseries have a tseries.name being cell
   for i = 1:num_tseries
      if iscell(data.tseries{i}.name)
         data.tseries{i}.name = data.tseries{i}.name{1};
      end
   end
   
   success = 1;
   tnow = datetime('now');
   str = sprintf( "\tSuccessfully opened voltage file %s at %s\n", fname, datestr(tnow) );
   printMessage('on',  'Keywords', str);
   
end

