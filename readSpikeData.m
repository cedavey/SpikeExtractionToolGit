% voltage = readSpikeData(filename, [type])
%  
% [voltage, time] = readSpikeData(filename)
%
% If type is set to something other than voltage then the SET GUI
% datastruct is returned instead of voltage. These files must be .mat or
% .smr or .txt
%
% Reads voltage & time data for darpa project from smr, txt, or mat file
function varargout = readSpikeData(filename, type)
   if nargin<2, type = []; end
   varargout = cell(1, nargout);
%     [varargout{:}] = deal([]); 
   data = [];
   % if filetype given extract type from filename
   type = [];
   if any('.'==filename)
      ind  = find(filename=='.', 1, 'last');
      type = ternaryOp(length(filename)>ind, filename(ind+1:end), []);
   end

    % if we know the filetype load it, else try different loading options &
    % see what works
    if ~isempty(type)
        switch type
            case 'smr'
                data = loadSMRfile(filename, true);
                if isempty(data)
                    str = sprintf('readSpikeData: could not load data, exiting...',filename);
                    cprintf('*Errors', str);
                    return;
                end
            case 'txt'
                data = loadTextSpike2file(filename, true);
            case 'mat'
                data = loadMatlabSpike2file(filename, true);
        end

    % try opening the file using diff file type options
    else
        try 
            data = loadMatlabSpike2file(filename, false);
        catch
            try
                 data = loadTextSpike2file(filename, false);
            catch
                try
                    data = loadSMRfile([filename '.smr'], false);
                catch 
                    str = sprintf('readSpikeData: file %s type is unrecognised, exiting...',filename);
                    cprintf('*Errors', str);
                    return;
                end
            end
        end
    end
    if isempty(data)
       varargout = deal({}); 
       return; 
    end
   
    % if data is a structure then this may be a previously generated SET
    % structure that has been saved - e.g. AP templates, spike rates etc   
    if isstruct(data)
       dname = fieldnames(data);
       SETdata  = data.(dname{1});
       % if we have a field 'type' then assume it's a SET struct
       if isfield(SETdata, 'type')
          if nargout>=1
             varargout{1} = SETdata;
          end
          if nargout>1
             varargout(2:nargout) = cell(length(2:nargout),1);
          end
          return;
       end
    end
    if nargout>=1
       if ~isempty( type )
          varargout{1} = data;
       else
          % try extracting as voltage but if it fails try extracting as data
          try
             varargout{1} = single(data.voltage); % prob already is single 
          catch
             varargout{1} = data.data;
          end
       end
    end
    if nargout>=2
        varargout{2} = single(data.time); % reduce storage cost
    end        

end

% Extract data from matlab file & look for time & vagus voltage
function data = loadTextSpike2file(filename, verbose)
    data  = [];
    fdata = tdfread(filename); % returns data in struct
    [time, voltage] = extractTimeVoltageFromStruct( fdata );
    
    data.time    = time; 
    data.voltage = voltage; 
    
    if verbose
        if isempty(time)
            str = 'Error finding time in data struct\n';
            cprintf('*Errors', str);
        end
        if isempty(voltage)
            str = 'Error finding time in data struct\n';
            cprintf('*Errors', str);
        end
    end
end


% Extract data from matlab file & look for time & vagus voltage
function data = loadMatlabSpike2file(filename, verbose)
   data = [];
   if exist(filename, 'file')
      try
         data   = load(filename); % returns data in struct
         if isempty(data), error('File is empty perhaps??'); end
      catch ME
         % assume we've run outta memory
         str = sprintf('\nI think you''ve run outta memory - gotta get smart\n\n');
         cprintf('Keywords*', str);
         return;
      end
   else
      str = sprintf('\n%s file not found, exiting...\n\n', filename);
      cprintf('Errors*', str); 
      return;
   end
   
   fnames = fields(data);
   % sometimes data is nestled one level deeper if saved as a struct
   if length( fnames ) == 1
      data = data.(fnames{1});
      fnames = fields(data);
   end
   if ~any( strcmpi( fnames, 'time' ) )
      if ~any( strcmpi( fnames, 'time' ) )      
         str = sprintf( '\n\tNo valid structure found, exiting\n');
         cprintf('Errors*', str);
         data = [];
      end
   end
   
end

function [data, ok] = loadSMRfile(filename, verbose)
   if nargin<2, verbose = true; end
   data = []; ok = false;
   if ~exist(filename, 'file')
      if verbose
         str = sprintf('%s file not found, exiting...\n', filename); 
         cprintf('*Errors', str);
      end
        return;
    end
   try
      fhand = CEDS64Open(filename);
   catch
      try
         % load CED path to manage with smr files
         cedpath = getenv('CEDS64ML'); % should hold the path to the CEDS64ML folder
         if isempty(cedpath)
            % if moved from default path, need to set the env variable accordingly
            % (this isn't setting to new path, just showing how to do it)
            cedpath = 'C:\CEDMATLAB\CEDS64ML'; % standard default CED path
            setenv('CEDS64ML', cedpath);       % do this if env var not set (correctly)
            if isempty(cedpath)
               if verbose
                  str = sprintf('Can''t find CED library, exiting...');
                  cprintf('*Errors', str);
                  return;
               end
            end
         end
         addpath( cedpath );
         CEDS64LoadLib( cedpath );
         fhand = CEDS64Open( filename );

      % unload the library when finished - I'ved moved this to closeGUI in
      % SpikeExtractionTool, just in case user's loading multiple files
      % unloadlibrary ceds64int;
      catch
         str = 'Can''t read smr files - need to load the (Windows only) CED library';
         cprintf('Error*', str);
         return;
      end
   end
   if (fhand <= 0)
      % unloadlibrary ceds64int; % moved to SpikExtractionTool->closeGui
      return; 
   end

   try
      % find channels available to write to
      maxchan    = CEDS64MaxChan( fhand );
      firstavail = CEDS64GetFreeChan( fhand );
         %   iType - 0 no channel (i.e. unused)
         %           1 Waveform channel
         %           2 Event (falling)
         %           3 Event (rising)
         %           4 Event (both)
         %           5 Marker
         %           6 Wavemark
         %           7 Realmark
         %           8 TextMark
         %           9 Realwave
      chantype = arrayfun( @(i) CEDS64ChanType( fhand, i ), firstavail:maxchan, 'uni', 1 );
      chantype = arrayfun( @(i) CEDS64ChanType( fhand, i ), 1:maxchan, 'uni', 1 );
      avail    = find( chantype~=0 );
      navail   = length( avail );

      list     = struct('type','list','format','integer','style', 'popupmenu','size',0);
      formats  = list;
      formats.limits = [1 navail];
      % numbers > 10 have 2 chars per number, else 1
      nchars   = ternaryOp( any( avail >= 10 ), 2, 1 );
      formats.items  = mat2cell( num2str( avail' ), ones(navail,1), nchars); % 2 chars per channel


      % need to update default value for lists from being string to
      % being index into list, else inputsdlg has a hissy
       def{1}  = 1;

      % ask user which channel they want to read from
      prompt   = {'Which channel would you like to read from?'};
      name     = 'Spike2 channel';
      % formats.size = -1 --> % auto-expand width and auto-set height
      options  = struct( 'Resize','on', 'WindowStyle','normal', 'Interpreter','tex' );
      default  = {1};   
      [answer,canceled] = inputsdlg(prompt, name, formats, default, options);
      if canceled
         [ ok ] = CEDS64Close( fhand );
         return;
      end
      channel = avail( answer{1} );

      
   catch
      [ ok ] = CEDS64Close( fhand );
   end
    
   % get waveform data from channel 1; +1 so the read gets the last data point 
   try
      maxTimeTicks = CEDS64ChanMaxTime( fhand, 1 ) + 1; 
      [numRead, voltage, fTime] = CEDS64ReadWaveF( fhand, channel, maxTimeTicks, 0, maxTimeTicks );
   catch ME
      str = sprintf( 'Error reading smr file %s on channel %d\n', filename, channel );
      cprintf('*Errors', str);
   end
   % returns negative value in fRead if error
   if verbose && numRead < 0 
      str = sprintf('Error reading smr file %s, returning empty data', filename);
      cprintf('*Errors', str);
   end
    
   % time resolution is time base of file * resolution of channel in ticks
   dt   = CEDS64ChanDiv(fhand, channel) *  CEDS64TimeBase(fhand);
   time = ( 1:numRead )' * dt;
   
   data.time  = time(:);
   data.data  = voltage(:);
   data.dt    = dt; 
   data.type  = 'voltage'; 
   data.params= []; 
   [ok,title] = CEDS64ChanTitle( fhand, channel ); % use chan title if it exists
   title      = ternaryOp( ~ok, 'voltage', title );
   data.name  = title; 
   
   [~, dataname, ext] = fileparts( filename ); 
   % for Martin! Because channel names are constant across files so he's
   % getting confused with multiple files open at once
   data.name  = dataname; 
    
   ok = CEDS64Close( fhand ); % if negative then error on closing
end


% Extract time & voltage from data struct by searching for field name with
% time & vagus in it, respectively
function [time, voltage] = extractTimeVoltageFromStruct( fdata )
    time    = extractDataFromStruct(data, 'time');
    voltage = extractDataFromStruct(data, 'vagus');
end

% Extract data from struct by identifying data name in struct fields
function x = extractDataFromStruct(data, dname)
    x      = [];
    if ~isstruct(data), return; end
    
    fnames = cellfun(@lower, fieldnames(fdata),'uni',0); % data fields in lower case
    ind    = cellfun(@any, strfind(fnames, dname));      % find vagus voltage field 
    
    % if can't find data name in struct fields then return with empty data
    if isempty(ind)
        return;
    end
    
    % if found more than 1 field in struct that matches data name, ask user
    if length(ind)>1
        str    = sprintf('More than 1 channel with %s in the name, select channel: \n', dname);
        inst   = arrayfun(@(i) sprintf('\tEnter %d for %s \n', i, fnames{ind(i)}),...
                          1:length(ind), 'uni', 0); 
        cprintf('Keywords', str);
        instr  = input(char(inst)');
        if isnumeric(instr) && (0<instr && instr<length(ind))
            ind = ind(instr);
        else
            str = sprintf('Invalid option, exiting with empty data...\n');
            cprintf('*Errors', str);
            return;
        end        
    end
    
    x = fdata.(fnames{ind});
end














