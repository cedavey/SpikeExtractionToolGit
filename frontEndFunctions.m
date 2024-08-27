classdef frontEndFunctions

   properties (Access = public)
      uiType % 'app' or 'gui'
   end

methods (Static)
   % Constructor
   function obj = frontEndFunctions(varargin)
      if nargin == 0, error('Not enough input arguments');end

      if nargin == 1
         obj.uiType = varargin{1};
      end
   end

   %% Shared methods
   function aboutMenuItem(h)
   % hObject    handle to aboutMenuItem (see GCBO)
   % eventdata  reserved - to be defined in a future version of MATLAB
   % handles    structure with handles and user data (see GUIDATA)
      aboutWindow = dialog();
      str = sprintf(['','\n',...
                     'The University of Melbourne','\n',...
                     '','\n',...
                     'Biomedical Engineering','\n',...
                     '','\n',...
                     'Dr. Katie Davey','\n',...
                     'Dr. Martin Stebbing','\n',...
                     'Dr. Artemio Soto Breceda','\n',...
                     '']);
      uicontrol('parent',aboutWindow,'Style','text',...
            'String', str,'Position',[187,0,187,200],'FontSize',11);

      % Get location of log files
      path = getFilePath();
      
      if strcmp('gui', h.f.uiType)
        logo = uicontrol('parent',aboutWindow,'Style','pushbutton',...
            'Position',[187,200,187,187]);
        [x,map]=imread([path 'fig' filesep 'unimelb.png']); % Load the zoom icon
        I2=imresize(x, [187 187]); % Resize icon
        logo.CData = I2; % Assign icon to the button
      else
         logo = axes('position',[0.3 0.475 0.375 0.475]);
         I = imread(['.' filesep 'fig' filesep 'unimelb.png']);
         imagesc(I);
         uistack(logo,'bottom');
         set(logo,'handlevisibility','off', ...
         'visible','off')
      end
   end

   function add_voltage(h)
      h.data.last_tseries = h.data.curr_tseries;
      new_num_tseries     = h.data.num_tseries + 1;

      % open image files & retrieve matrices
      [h.data, success]   = openVoltageFile(h.data);

      h = h.f.resetSliders(h);

      if success~=1
         str = sprintf('Error opening file, ignoring');
         displayErrorMsg(str);
         return
      end
      % loading new voltage tseries succeeded

      if strcmp('gui', h.f.uiType)
         % udpate drop down lists - gotta be done for success = 0 or 1 (?)
         set(h.curr_signal, 'String', h.data.tseries_str);
         % set voltage to 1st newly loaded voltage tseries
         set(h.curr_signal, 'Value', new_num_tseries);
      else
         h.curr_signal.Items = h.data.tseries_str;
         set(h.curr_signal, 'Value', h.curr_signal.Items(new_num_tseries));
      end
      h.f.curr_signal(h.curr_signal, 'add_voltage', h);
   end

   function h = automatic_params(h, hObject)
      if hObject.Value
         h.set_tool_params.Enable = 'off';
         h.options.auto_params    = true;
      else
         h.set_tool_params.Enable = 'on';
         h.options.auto_params    = false;
      end
      if strcmp('gui', h.f.uiType), guidata(hObject,h); end% saves changes to handles
   end

   function clear_voltages(hObject, h, vtc)
      try
         vtc = vtc.Value;
         delete(hObject.Parent);
         for i = numel(vtc):-1:1
            h = h.f.removeVoltage(h, vtc(i));
         end
      catch E
         caught_error = runtimeErrorHandler(E);
         if ~isempty(caught_error), rethrow(caught_error); end
      end
   end
   
   function runBatchProcessing(handles, varargin)
      if nargin > 1
         hObject = varargin{1};
         eventdata = varargin{2};
      end
      
      try
         opts = batchProcessing();
      catch E
         if strcmp('No_choice', E.message)
            str = sprintf('\tBatch process cancelled\n');
            printMessage('off', 'Error', str);
            return;
         else
            caught_error = runtimeErrorHandler(E, 'rethrow');
            if ~isempty(caught_error), rethrow(caught_error); end
         end
      end

      % Load files
      for i = 1:numel(opts.files)
         last_dir = handles.data.last_dir; % keep before we write over it
         old_numtseries = handles.data.num_tseries;
         % get rid of all previous data
         handles = toggleSETGUIstate(handles,'off');
         handles.data.last_dir = last_dir;
         data    = handles.data;    % get user data from gui handle
         [data, success] = openVoltageFile(data, 'batch', opts.path, opts.files(i));
         data.last_tseries = 1;
         data.curr_tseries = 1;
         data.last_tool    = 1;
         data.curr_tool    = 1;

         if success==0 % if success==0 --> no valid images found or user cancelled
            displayErrorMsg('No valid voltage data found - please reload');
            % don't update handles with the data_struct changes
            return;
         elseif success==-1
            if old_numtseries>0 % user cancelled out of open file dialogue
               displayErrorMsg('Load voltage cancelled, but old data was removed (sorry!)');
            end
            return;
         end
         % if new data loaded re-enable GUI
         handles.data = data;
         handles = toggleSETGUIstate(handles,'on');
                  
         if strcmp('gui',handles.f.uiType)
            guidata(hObject,handles);  % saves the change to handles
            set(handles.curr_signal, 'String', data.tseries_str);
            set(handles.curr_signal, 'Value',  1);
         else
            set(handles.curr_signal, 'Items', data.tseries_str);
            set(handles.curr_signal, 'Value',  handles.curr_signal.Items(1));
         end

         % curr_signal_changed doesn't return handles so we have to save handles
         % manually, then request a fresh copy using guidata
         % curr_signal_changed(handles.curr_signal, '', handles);
         handles.f.curr_signal(handles.curr_signal, '', handles);

         [tl, ml] = getBatchToolList(); % Load the list of tools available
         % Tool method and params
         for j = 1:numel(opts.tool)
            % Implement the tool using the method & params requested
            if strcmp('gui',handles.f.uiType)
               set(handles.tool_list, 'String', opts.tool(j));
               set(handles.tool_list, 'Value', 1);
            else
               set(handles.tool_list, 'Items', opts.tool);
               set(handles.tool_list, 'Value',  handles.tool_list.Items(1));
            end
            idx = find(ismember(tl, opts.tool{j}));
            if ~isempty(idx)
               method = ml{idx};
            else
               error('Wrong tool selected for batch processing.');
            end
            
            if strcmp('gui',handles.f.uiType)
               set(handles.method_list, 'String',{method});
               set(handles.method_list, 'Value', 1);
            else
               set(handles.method_list, 'Items', {method});
               set(handles.method_list, 'Value',  handles.method_list.Items(1));
            end
            % Automatic parameters will always be true for batch processing
            handles.options.auto_params = true;
            handles.options.isBatch     = true;
            handles.options.batchPath   = opts.saveFolder;
            handles.options.debugOption = 'none';

            try
               handles = handles.f.run_tool(handles);
               % Save all open tseries. These will be all the processed tools on the
               % current file
               handles.f.save_voltage(handles, hObject);
               handles.options.isBatch = false;
            catch E
               str = sprintf('\tBatch processing failed at i = %s, j = %s\n\tError: %s\n',opts.files{i}, opts.tool{j}, E.message);
               caught_error = runtimeErrorHandler(E, 'message', str);
               if ~isempty(caught_error), rethrow(caught_error); end
               handles.options.isBatch = false;
               if strcmp('gui', handles.f.uiType), guidata(hObject,handles); end% saves changes to handles
               break;
            end
         end % Tools and Params <j>
      end % Files <i>
      if strcmp('gui', handles.f.uiType), guidata(hObject,handles); end% saves changes to handles
   end

   function curr_signal(hObject, eventdata, handles)
      if ~handles.f.haveUserData(handles)
         return;
      end
      % find out what the user has selected to view
      [tseries, ts_num, data_type] = handles.f.getCurrentVoltage(handles);
      % add current tseries to last viewed tseries before it's overwritten
      if ~strcmp('add_voltage', eventdata)
         % Only if this is not a new voltage
         handles.data.last_tseries    = handles.data.curr_tseries;
      end
      % update current tseries to tseries chosen by user
      handles.data.curr_tseries    = ts_num;

      % if data type of tseries hasn't changed then tool list doesn't need to
      % change, but if we're displaying a new data type then the tool list and
      % available methods for the tool needs to change
      last_tseries = handles.data.tseries{handles.data.last_tseries};
      last_type    = handles.data.tseries{handles.data.last_tseries}.type;
      [last_tlim, last_vlim] = getTimeAndVoltageLimits(last_tseries);
      [curr_tlim, curr_vlim] = getTimeAndVoltageLimits(tseries);

      if strcmp('load_voltage', eventdata) || ~strcmpi(data_type, last_type)
         % only remember last tool for voltage data cuz not many options for the others
         tool_num    = ternaryOp( strcmpi(data_type, 'voltage'), handles.data.last_tool, 1);
         if (isempty(tool_num) || tool_num==0), tool_num=1; end
         tool_list   = getSETToolList(data_type);
         tool        = tool_list{tool_num};
         method_list = getSETToolMethodsList(tool, data_type);
         if strcmp('gui', handles.f.uiType)
            set(handles.tool_list,   'String', tool_list);
            set(handles.tool_list,   'Value',  tool_num);
            set(handles.method_list, 'String', method_list);
            set(handles.method_list, 'Value', 1);
         else
            set(handles.tool_list,   'Items', tool_list);
            set(handles.tool_list,   'Value',  tool_list(tool_num));
            set(handles.method_list, 'Items', method_list);
            set(handles.method_list, 'Value', method_list(1));
         end
      end

      % update voltage & time sliders - keep them the same if timeseries length
      % is the same as the last viewed tseries
      if ~compareFloats(curr_tlim(2), last_tlim(2), 0.01, 'perc')
         % set time slider and time max/min text boxes from time length of data
         handles.data.tlim = curr_tlim;
         mint = handles.data.tlim(1); maxt = handles.data.tlim(2);
         set(handles.time_slider, 'Value',  0);
         if strcmp('gui', handles.f.uiType)
            set(handles.time_max,    'String', sprintf('%.2f',maxt));
            set(handles.time_min,    'String', sprintf('%.2f',mint));
         else
            set(handles.time_max,    'Value', sprintf('%.2f',maxt));
            set(handles.time_min,    'Value', sprintf('%.2f',mint));
         end
      end

      if strcmp('load_voltage', eventdata) || ~compareFloats(curr_vlim(2), last_vlim(2), 0.01, 'perc') || ts_num ~= handles.data.last_tseries
         % set time slider and time max/min text boxes from scale of data
         handles.data.vlim = curr_vlim;
         minv = handles.data.vlim(1); maxv = handles.data.vlim(2);

         if strcmp('gui', handles.f.uiType)
            set(handles.voltage_max,    'String', sprintf('%.2f',maxv));
            set(handles.voltage_min,    'String', sprintf('%.2f',minv));

            % udpate figure to show tseries chosen by user
            guidata(hObject,handles); % saves changes to handles
            handles = updateSETFigure(handles, tseries);
            updateGUIParams(handles, tseries);

            guidata(hObject,handles); % saves changes to handles
         else
            set(handles.voltage_max,    'Value', sprintf('%.2f',maxv));
            set(handles.voltage_min,    'Value', sprintf('%.2f',minv));

            % udpate figure to show tseries chosen by user
            guidata(hObject,handles); % saves changes to handles
            handles = updateSETFigure(handles, tseries);
            updateGUIParams(handles, tseries);
         end
      end
   end

   function [tseries, ts_num, type, ts_name] = getCurrentVoltage(handles)
   % extract time series data, index number into tseries cell array, data type
   % of time series, and the name of the time series
   % if try fails then we're doing if for viewVolume_newWindow - get stat
   % from data struct instead
   try
      if strcmp('gui', handles.f.uiType)
         ts_num  = get(handles.curr_signal,'Value');
      else
         v = cellfun(@(x) strcmp(x,handles.curr_signal.Value), handles.curr_signal.Items);
         ts_num = find(v == 1);
      end
      tseries = handles.data.tseries{ts_num};
      ts_name = handles.data.tseries_str{ts_num};
      if iscell(ts_name)
        ts_name = ts_name{1};
      end
      if iscell(tseries.name)
        tseries.name = tseries.name{1};
      end
      type    = tseries.type;
   catch ME
      ts_num  = handles.data.curr_tseries;
      tseries = handles.data.tseries{ts_num};
      ts_name = handles.data.tseries_str{ts_num};
      if iscell(ts_name)
        ts_name = ts_name{1};
      end
      if iscell(tseries.name)
        tseries.name = tseries.name{1};
      end
      type    = tseries.type;
   end
   end

   function name = getFileName(instruct, defVal, maxLength)
   % prompt user for name of volume saved to file (so they're not stupidly
      % long)
      if isempty(defVal)
         defVal = 'image';
      end
      name = inputdlg(instruct, 'Variable name', 1,{defVal});
      if isempty(name) % user's cancelled process
         return
      end
      name = name{1};
      if exist('maxLength','var') && ~isempty(maxLength)
         while length(name)>maxLength
            displayErrorMsg(['Can''t exceed ' num2str(maxLength) ' chars, truncating']);
            name = inputdlg(instruct(1:maxLength), 'Variable name', 1,{defVal});
         end
      end
   end

   function userids = getUserTemplateDeleteIDs(tseries, method_params)
      % user configuration parameters:
      %   ap.merge_templates.user_selection.template_to_merge_with.value = 1;
      %   ap.merge_templates.user_selection.number_of_templates_to_merge.value   = 1;
      userids = [];
      nap     = size( tseries.data, 2 );
      ndelete  = method_params.number_of_templates_to_remove.value;
      % make sure there's enough templates to merge as many as user requested
      if ndelete >= nap
         str = 'You''re trying to delete more templates than you actually have, try again';
         displayErrorMsg(str);
         return;
      end

      % get IDs of templates available to merge
      ids = 1:nap;
      list = struct('type','list', 'format','integer', 'style','popupmenu', 'size',0);
      for ii=1:ndelete
         def{ii,1}    = ids( ii );
         descript     = sprintf( 'Delete %d templates');
         name         = sprintf( 'ap%d', ii ); % ids( ii ) );
         units        = 'integer';
         prompt{ii,1} = sprintf( 'Select template ID to remove');
         prompt{ii,2} = name;    % name of struct field for result
         prompt{ii,3} = units;   % units of parameter

         % list options
         options      = ids;
         tmp          = list;
         tmp.limits   = [1 length(options)];
         tmp.items    = options;
         formats(ii,1)= tmp;

         % need to update default value for lists from being string to
         % being index into list, else inputsdlg has a hissy
         def{ii,1}    = find( def{ii,1} == options );
      end
      other.Resize      = 'on';
      other.WindowStyle = 'normal';
      other.Interpreter = 'tex';

      dlg_title = sprintf( 'Select AP templates to remove');
      try
         [userparams, cancel] = inputsdlg(prompt, dlg_title, formats, def, other);
      catch ME
         str = getCatchMEstring( ME, 'Error setting parameters', false );
         displayErrorMsg( 'Error setting parameters, reverting to old values' );
         userids = [];
         caught_error = runtimeErrorHandler(ME,'ignore');
         if ~isempty(caught_error), rethrow(caught_error); end
         return;
      end
      if cancel
         userids = [];
         return;
      end

      % replace list indices with template ids
      names = fieldnames (userparams );
      for i=1:ndelete
         % for the list types extract value from index into list
         userparams.( names{i} ) = formats(i).items( userparams.( names{i} ) );
      end

      % if any template ids are repeated, delete repeat
      userids = cellfun( @(f) userparams.(f), names );
      if length( unique( userids ) ) ~= length( userids )
         displayErrorMsg( 'Ignoring repeat template IDs' );
      end
      userids = unique( userids );
   end

   function userids = getUserTemplateMergeIDs(tseries, method_params)
      % user configuration parameters:
      %   ap.merge_templates.user_selection.template_to_merge_with.value = 1;
      %   ap.merge_templates.user_selection.number_of_templates_to_merge.value   = 1;
      userids = [];
      nap     = size( tseries.data, 2 );
      mergeid = method_params.template_to_merge_with.value;
      nmerge  = method_params.number_of_templates_to_merge.value;
      % make sure there's enough templates to merge as many as user requested
      if nmerge >= nap
         str = 'You''re trying to merge more templates than you actually have, try again';
         displayErrorMsg(str);
         return;
      end

      % get IDs of templates available to merge
      ids = 1:nap;
      ids( mergeid ) = []; % remove template we're merging into from list
      list = struct('type','list', 'format','integer', 'style','popupmenu', 'size',0);
      for ii=1:nmerge
         def{ii,1}    = ids( ii );
         descript     = sprintf( 'Template number to merge in to template %d', mergeid );
         name         = sprintf( 'ap%d', ii ); % ids( ii ) );
         units        = 'integer';
         prompt{ii,1} = sprintf( 'Select template ID of AP to merge into template %d', mergeid );
         prompt{ii,2} = name;    % name of struct field for result
         prompt{ii,3} = units;   % units of parameter

         % list options
         options      = ids;
         tmp          = list;
         tmp.limits   = [1 length(options)];
         tmp.items    = options;
         formats(ii,1)= tmp;

         % need to update default value for lists from being string to
         % being index into list, else inputsdlg has a hissy
         def{ii,1}    = find( def{ii,1} == options );
      end
      other.Resize      = 'on';
      other.WindowStyle = 'normal';
      other.Interpreter = 'tex';

      dlg_title = sprintf( 'Select AP templates to merge into template %d', mergeid );
      try
         [userparams, cancel] = inputsdlg(prompt, dlg_title, formats, def, other);
      catch ME
         str = getCatchMEstring( ME, 'Error setting parameters', false );
         displayErrorMsg( 'Error setting parameters, reverting to old values' );
         userids = [];
         caught_error = runtimeErrorHandler(ME,'ignore');
         if ~isempty(caught_error), rethrow(caught_error); end
         return;
      end
      if cancel
         userids = [];
         return;
      end

      % replace list indices with template ids
      names = fieldnames (userparams );
      for i=1:nmerge
         % for the list types extract value from index into list
         userparams.( names{i} ) = formats(i).items( userparams.( names{i} ) );
      end

      % if any template ids are repeated, delete repeat
      userids = cellfun( @(f) userparams.(f), names );
      if length( unique( userids ) ) ~= length( userids )
         displayErrorMsg( 'Ignoring repeat template IDs' );
      end
      userids = unique( userids );
   end

   function have = haveUserData(handles)
      have = 0;
      % if no timeseries to display don't try moving crosshairs etc
      try
         if handles.data.num_tseries==0
            return;
         end
      catch   % no data_struct (error trying to access it)
         return;
      end
      have = 1;
   end

   function load_voltage(varargin)
      if nargin == 1
         h = varargin{1};
      else
         hObject = varargin{1};
         h = varargin{2};
      end

      last_dir = h.data.last_dir; % keep before we write over it
      old_numtseries = h.data.num_tseries;

      % get rid of all previous data
      h = toggleSETGUIstate(h,'off');
      h.data.last_dir = last_dir; % add last dir back to handle

      if strcmp('gui',h.f.uiType)
         guidata(hObject,h);  % saves the change to handles
      end
      data    = h.data;    % get user data from gui handle

      % open smr, txt, or mat file

      try
        temp_dir = data.last_dir;
        data.last_dir = userpath;
        [data, success] = openVoltageFile(data);  
      catch 
        data.last_dir = temp_dir;
        [data, success] = openVoltageFile(data);  
      end

      data.last_tseries = 1;
      data.curr_tseries = 1;
      data.last_tool    = 1;
      data.curr_tool    = 1;

      if success==0 % if success==0 --> no valid images found or user cancelled
         displayErrorMsg('No valid voltage data found.');
         % don't update handles with the data_struct changes
         return;
      elseif success==-1
         if old_numtseries>0 % user cancelled out of open file dialogue
            displayErrorMsg('Load voltage cancelled, but old data was removed (sorry!)');
         end
         return;
      end
      % if new data loaded re-enable GUI
      h.data = data;
      h = toggleSETGUIstate(h,'on');

      h = h.f.resetSliders(h);

      if strcmp('gui',h.f.uiType)
         guidata(hObject,h);  % saves the change to handles
         set(h.curr_signal, 'String', data.tseries_str);
         set(h.curr_signal, 'Value',  1);
         guidata(hObject,h);
      else
         set(h.curr_signal, 'Items', data.tseries_str);
         set(h.curr_signal, 'Value',  h.curr_signal.Items(1));
      end

      % curr_signal_Callback doesn't return handles so we have to save handles
      % manually, then request a fresh copy using guidata
      h.f.curr_signal(h.curr_signal, 'load_voltage', h);
      if strcmp('gui',h.f.uiType), guidata(hObject,h); end
   end

   % create a new gui using the current data displayed
   function new_figure(h)
      new_gui     = SpikeExtractionTool;
      new_handles = guidata(gcf);
      tseries     = h.f.getCurrentVoltage(h);
      h.data.guihandles = [h.data.guihandles; new_gui.guihandles]; % record new gui in our data struct
      guidata(h.figure1, handles);

      %% Generate user data object to store
      new_gui.num_tseries = 1;
      new_gui.tseries     = {tseries};             % initialise cell array of time series
      new_gui.curr_tseries= 1;                     % current time series
      new_gui.tseries_str = {tseries.name};        % initialise names of time series
      new_gui.used_names  = {{tseries.name}, 1};   % list of names already used by user
      new_gui.last_tseries= 1;   % last time series - restore if current ts deleted
      new_gui.tlim        = h.data.tlim;
      new_gui.vlim        = h.data.vlim;
      new_gui.last_tool   = 1;                     % record of user's last tool
      new_gui.params      = getDefaultToolParams;  % default params for all tools & implementation methods
      new_gui.last_dir    = h.data.last_dir;       % where they opened gui from
      % copy of all SET gui handles - for each gui, make sure it's own handle is 1st
      new_gui.guihandles  = [new_gui.guihandles; h.data.guihandles(1:end-1)];
      new_app.data        = new_gui;               % record user data in handle
      new_handles         = toggleSETGUIstate(new_handles, 'on'); % switch everything off until data's loaded

      % udpate voltage timeseries drop down list
      set(new_app.curr_signal, 'String', new_app.data.tseries_str);
      set(new_app.curr_signal, 'Value', 1);

      % plot data that was put in new figure
      % curr_signal_Callback( new_app.curr_signal, [], new_handles );
      new_app.f.curr_signal(new_app.curr_signal, [], new_handles);
   end

   function varargout = removeVoltage(h, varargin)
      if (numel(varargin) == 0) || cellfun(@isempty, varargin)
         remove_tseries = h.data.curr_tseries;
         varargin = [];
         varargout = {};
      else
         remove_tseries = varargin{1}; % In case we are removing a voltage different to the currently active
      end
         tseries_str    = h.data.tseries_str{remove_tseries};
         last_tseries   = h.data.last_tseries;
         num_tseries    = h.data.num_tseries;

      % if the last remaining timeseries is being deleted, clear gui
      if num_tseries==1
         h = toggleSETGUIstate(h, 'off');
         guidata(h.figure1, h);
         if numel(varargin) == 1,varargout = {h};end % When removing multiple voltages at once, we want to update the state of handles for the caller method
         return;

         % last time series may be the same as current time series (being removed)
         %  after deleting a timeseries or some such thing, so show something else
      elseif h.data.last_tseries==remove_tseries
         indices     = 1:num_tseries;
         indices(remove_tseries) = [];
         h.data.last_tseries = indices(1);
      end
      % display last tseries before removing this one so we can check if data
      % types are the same etc
      if strcmp('gui', h.f.uiType)
         set(h.curr_signal, 'Value', h.data.last_tseries);
         % Update figure & handles structure
         h.f.curr_signal(h.curr_signal, [], h);
         h = guidata(h.figure1);
      else
         set(h.curr_signal, 'Value', h.curr_signal.Items(h.data.last_tseries));
         % Update figure & handles structure
         h.f.curr_signal(h.curr_signal, [], h);
      end
      % curr_signal changes last & current timeseries indices, so change back
      % (we wanted the function to display the previous timeseries, but the
      % user didn't actually select it themselves so don't update current &
      % last timeseries indices)
      h.data.last_tseries = last_tseries;
      h.data.curr_tseries = remove_tseries;

      % now do the actual removing
      h.data.used_names = removeStringFromList(h.data.used_names,  tseries_str);
      h.data.tseries(remove_tseries)     = [];
      h.data.tseries_str(remove_tseries) = [];
      h.data.num_tseries              = h.data.num_tseries - 1;

      if h.data.num_tseries==0
         h = toggleGUIstate(h,'off');
         if numel(varargin) == 1,varargout = {h};end % When removing multiple voltages at once, we want to update the state of handles for the caller method
         return
      end
      % update our indices for last tseries if necessary
      if h.data.last_tseries == remove_tseries
         h.data.last_tseries = 1;
      elseif h.data.last_tseries > remove_tseries
         h.data.last_tseries = h.data.last_tseries - 1;
      end
      if h.data.last_tseries > length(h.data.tseries_str)
         str = sprintf('last timeseries index (%d) larger than number of strings (%d)',...
            h.data.last_tseries, length(h.data.tseries_str) );
         displayErrorMsg(str);
         h.data.last_tseries = 1;
      end
      h.data.curr_tseries = h.data.last_tseries;
      if strcmp('gui', h.f.uiType)
         set(h.curr_signal, 'String', h.data.tseries_str);
         set(h.curr_signal, 'Value', h.data.last_tseries);
         h.data.curr_tseries = h.data.last_tseries;
         guidata(h.figure1, h);
      else
         set(h.curr_signal, 'Items', h.data.tseries_str);
         set(h.curr_signal, 'Value', h.curr_signal.Items(h.data.last_tseries));
         h.data.curr_tseries = h.data.last_tseries;
      end

      if numel(varargin) == 1,varargout = {h};end % When removing multiple voltages at once, we want to update the state of handles for the caller method
   end
   
   function resetButtonClicked(h, hObject)
      tseries = h.data.tseries{h.data.curr_tseries};
      [curr_tlim, curr_vlim] = getTimeAndVoltageLimits(tseries);
      
      % set time slider and time max/min text boxes from time length of data
      h.data.tlim = curr_tlim;
      mint = h.data.tlim(1); maxt = h.data.tlim(2);
      set(h.time_slider, 'Value',  0);
      if strcmp('gui', h.f.uiType)
         set(h.time_max,    'String', sprintf('%.2f',maxt));
         set(h.time_min,    'String', sprintf('%.2f',mint));
      else
         set(h.time_max,    'Value', sprintf('%.2f',maxt));
         set(h.time_min,    'Value', sprintf('%.2f',mint));
      end
      
     % set time slider and time max/min text boxes from scale of data
      h.data.vlim = curr_vlim;
      minv = h.data.vlim(1); maxv = h.data.vlim(2);

      if strcmp('gui', h.f.uiType)
         set(h.voltage_max,    'String', sprintf('%.2f',maxv));
         set(h.voltage_min,    'String', sprintf('%.2f',minv));

         % udpate figure to show tseries chosen by user
         guidata(hObject,h); % saves changes to handles
         h = updateSETFigure(h, tseries);
         updateGUIParams(h, tseries);

         guidata(hObject,h); % saves changes to handles
      else
         set(h.voltage_max,    'Value', sprintf('%.2f',maxv));
         set(h.voltage_min,    'Value', sprintf('%.2f',minv));

         % udpate figure to show tseries chosen by user
         guidata(hObject,h); % saves changes to handles
         h = updateSETFigure(h, tseries);
         updateGUIParams(h, tseries);
      end
      
      h.f.resetSliders(h, hObject);
   end

   function h = resetSliders(h, varargin)
      if strcmp('gui', h.f.uiType) && (nargin > 1)
        hObject = varargin{1};
      end
      
      h.data.time_slider_max = false;

      h.data.zoomPercentage = [0 0]; % Records currently chosen zoom value
      h.data.displacementPercentage = [0 0.5]; % Records currently chosen displacement value
      % Get location of GUI files
      path = getFilePath();

      [x,~]=imread([path 'fig' filesep 'magnifierIcon.png']);% Load the zoom icon
      if strcmp('gui', h.f.uiType)
        I2=imresize(x, [28.6 85.8]); % Resize icon
        h.toggleZoomButton.CData = I2; % Assign icon to the button
        h.toggleZoomButton.UserData = 'zoom'; % Change state to zoom
        set(h.zoom_out_label,'Visible', 'on');
        set(h.zoom_in_label,'Visible', 'on');
      else
        h.toggleZoomButton.Icon = [path 'fig' filesep 'magnifierIcon.png'];
      end

      h.time_slider.Value = h.data.zoomPercentage(1); % Update the position of the slider to represent zoom.
      h.voltage_slider.Value = h.data.zoomPercentage(2); % Update the position of the slider to represent zoom.

      if strcmp('gui', h.f.uiType)
        h.time_slider.SliderStep    = [0.01 0.1]; % Make sure it is within 0 and 1.
        h.voltage_slider.SliderStep = [0.01 0.1];
      end

      if exist('hObject')
        % Update tooltip
        h = setTooltips(h, {hObject.Tag}, getTooltips({'reset_button'}));
      elseif isfield(h, 'reset_button')
        h = setTooltips(h, {h.reset_button.Tag}, getTooltips({'reset_button'}));
      end
      
      if strcmp('gui', h.f.uiType)
        guidata(h.reset_button, h);% guidata(hObject,h); 
      end

      h.f.time_slider_updated(h, h.time_slider);
      
   end

   function varargout = run_tool(h)
   printcol = 'comments*'; % message category dictates printing colour in cprintf
     tic;
     varargout = {};
      % Implement the tool using the method & params requested
      if strcmp('gui', h.f.uiType)
         tool_list     = get(h.tool_list, 'String');
         tool_num      = get(h.tool_list, 'Value');
         tool          = tool_list{tool_num};

         method_list   = get(h.method_list, 'String');
         method_num    = get(h.method_list, 'Value');
         method        = method_list{method_num};
      else
         tool_list     = get(h.tool_list, 'Items');
         tn            = cellfun(@(x) strcmp(x,h.tool_list.Value), h.tool_list.Items);
         tool_num      = find(tn == 1);
         tool          = tool_list{tool_num};

         method_list   = get(h.method_list, 'Items');
         tn            = cellfun(@(x) strcmp(x,h.method_list.Value), h.method_list.Items);
         method_num      = find(tn == 1);
         method        = method_list{method_num};
      end

      [tseries, ~, type] = h.f.getCurrentVoltage(h);
      if ~strcmpi('export to excel', tool)
        % Because 'export to excel' does not have params, only do this for the
        % other tools
        method_params = getToolAndMethodParams(h.data.params, type, tool, method);
      end

      switch lower(type)
         case 'voltage'
            % only remember last tool for voltage coz others don't have many options
            h.data.last_tool = h.data.curr_tool;
            if strcmp('gui', h.f.uiType)
               h.data.curr_tool = get(h.select_tool, 'Value');
            else
               ct = cellfun(@(x) strcmp(x,h.tool_list.Value), h.tool_list.Items);
               h.data.curr_tool = find(ct == 1);
            end
            switch lower(tool)
               case 'rescale'
                  str = sprintf( 'Rescaling voltage\n' );
                  cprintf( printcol, str );
                  if ~strcmp('separate',method_params.select_peaks.value)
                     [voltage, Rest] = rescaleVoltage(tseries, method, method_params, h.options);
                  else
                     method_params.select_peaks.value = 'positive';
                     [voltage, Rest, prev_params] = rescaleVoltage(tseries, method, method_params, h.options);
                     voltage = voltage.*heaviside(voltage);
                     method_params.select_peaks.value = 'negative';
                     [nvoltage, Nest] = rescaleVoltage(tseries, method, method_params, h.options);
                     nvoltage = -nvoltage;
                     nvoltage = nvoltage.*heaviside(nvoltage);
                     voltage = voltage - nvoltage;
                  end
                  if isempty(voltage)
                     return;
                  end
                  new_tseries.type    = 'voltage';
                  new_tseries.data    = voltage;
                  new_tseries.time    = tseries.time;
                  new_tseries.dt      = tseries.dt;
                  new_tseries.params  = method_params;
                  new_tseries.Rest    = Rest;
                  new_tseries.params.tool   = tool;
                  new_tseries.params.method = method;
                  tool_str            = [tseries.name '_' tool];
                  instruct            = ['Creating ' tool_str ': rename?'];


               case 'denoise'
                  str = sprintf( 'Denoising voltage\n' );
                  cprintf( printcol, str );
                  voltage = denoiseVoltage(tseries, method, method_params);
                  if isempty(voltage)
                     return;
                  end
                  new_tseries.type    ='voltage';
                  new_tseries.data    = voltage;
                  new_tseries.time    = tseries.time;
                  new_tseries.dt      = tseries.dt;
                  new_tseries.params  = method_params;
                  new_tseries.params.tool   = tool;
                  new_tseries.params.method = method;
                  tool_str            = [tseries.name '_' tool];
                  instruct            = ['Creating ' tool_str ': rename?'];

               case 'identify ap templates'
                  str = sprintf( 'Identifying action potential templates from timeseries\n' );
                  cprintf( printcol, str );
                  % get time series to apply AP templates to
                  [APtemplates, APfamily,alignment_inds_templates] = identifyTemplates(tseries, method, method_params,  h.options);
                  if isempty(APtemplates)
                     str = sprintf('\tNo templates could be identified.\n');
                     printMessage('off', 'Error', str);
                     return;
                  end
                  new_tseries.type    = 'ap';
                  new_tseries.data    = APtemplates;
                  % coded to allows different length spiking templates
                  if iscell( APtemplates )
                     for ti=1:length(APtemplates)
                        new_tseries.time{ti} = (1:size(APtemplates{ti},1))'*tseries.dt;
                     end
                  else
                     new_tseries.time    = (1:size(APtemplates,1))'*tseries.dt;
                  end
                  new_tseries.dt      = tseries.dt;
                  new_tseries.params  = method_params;
                  new_tseries.params.tool   = tool;
                  new_tseries.params.method = method;
                  new_tseries.APfamily= APfamily;
                  new_tseries.alignment_inds = alignment_inds_templates;
                 
                  tool_str            = [tseries.name '_APs'];
                  instruct            = ['Creating ' tool_str ': rename?'];

               case 'extract spikes'
                  % allows user to generate spikes directly from voltage
                  % timeseries, by generating AP templates enroute
                  if strcmpi(method,'k means')
                     str = sprintf( 'Extracting spikes from timeseries using k-means\n' );
                     cprintf( printcol, str );
                     % Extracting spikes directly from voltage through K-means
                     % clustering. Not using APs
                     [APspikes, APtimes, APfamily] = extractSpikesUsingKmeans(tseries, method_params, h.options.debugOption);
                     if isempty(APspikes)
                        return;
                     end
                     new_tseries.type    = 'spike';
                     new_tseries.data    = APspikes;
                     new_tseries.time    = tseries.time;
                     new_tseries.dt      = tseries.dt;
                     new_tseries.params  = method_params;
                     new_tseries.params.tool   = tool;
                     new_tseries.params.method = method;
                     new_tseries.APfamily= APfamily; % Here the AP templates are the same as the different axons
                     new_tseries.APstimes= APtimes; % record spike times for getting rate later
                     tool_str            = [tseries.name '_spikes'];
                     instruct            = ['Creating ' tool_str ': rename?'];

                  else
                     str = sprintf( 'Extracting spikes from timeseries by identifying action potential templates\n' );
                     cprintf( printcol, str );
                     [APtemplates, APfamily]= identifyTemplates(tseries, 'threshold', method_params);
                     if isempty(APtemplates)
                        str = sprintf('No templates found\n');
                        return;
                     end
                     APnumsamples        = cellfun(@(f) size(f,2), APfamily);
                     normAPs             = method_params.normalise_aps.value; % separate for when gen APs separately
                     [APspikes,APtimes]  = extractSpikesUsingTemplates( APtemplates, APnumsamples, tseries, method, method_params, normAPs, h.options );
                     if isempty(APspikes)
                        str = sprintf('\tNo axon families could be identified.\n');
                        printMessage('off', 'Error', str);
                        return;
                     end
                     new_tseries.type    = 'spike';
                     new_tseries.data    = APspikes;
                     new_tseries.time    = tseries.time;
                     new_tseries.dt      = tseries.dt;
                     new_tseries.params  = method_params;
                     new_tseries.params.tool   = tool;
                     new_tseries.params.method = method;
                     new_tseries.APfamily= APfamily;
                     new_tseries.APstimes= APtimes; % record spike times for getting rate later
                     tool_str            = [tseries.name '_spikes'];
                     instruct            = ['Creating ' tool_str ': rename?'];
                  end

               case 'utilities'
                  str = sprintf( 'Applying utilities\n' );
                  cprintf( printcol, str );
                  new_tseries = voltageUtilities(tseries, method, method_params);
                  if isempty(tseries)
                     return;
                  end
                  new_tseries.type    ='voltage';
                  new_tseries.name    = tseries.name;
                  new_tseries.params  = method_params;
                  new_tseries.params.tool   = tool;
                  new_tseries.params.method = method;
                  tool_str            = [tseries.name '_' method];
                  instruct            = ['Creating ' tool_str ': rename?'];

               otherwise
                  displayErrorMsg('This tool don''t exist, give up');
                  return;
            end

            %% Action potentials
         case 'ap'
            switch lower(tool)
               case 'extract spikes'
                  % if user hasn't set parameters explicitly they wouldn't
                  % have chosen a voltage timeseries to apply the AP templates to                  
                  if ~isfield(method_params, 'voltage_timeseries')
                     if ~isfield(h.options, 'isBatch') || ~h.options.isBatch 
                        % if user hasn't selected which voltage timeseries to
                        % extract spikes from, if there's only 1 timeseries
                        % there assume it's that one, else ask them to choose
                        % Only do this if Batch is not selected
                        data_types = getStructFieldFromCell( h.data.tseries, 'type' );
                        if length( data_types ) == 2 && any( strcmpi( 'voltage', data_types ) )
                           vind = find( strcmpi( 'voltage', data_types ) );
                           method_params.voltage_timeseries = struct;
                           method_params.voltage_timeseries.value = h.data.tseries{vind}.name;
                           method_params.voltage_timeseries.name = 'voltage timeseries';
                           method_params.voltage_timeseries.type = 'string';
                           
                        else
                           % if only 1 voltage data type use it for extracting spikes
                           dtypes = getStructFieldFromCell( h.data.tseries, 'type' );
                           vdata = strcmpi( 'voltage', dtypes );
                           if ( sum( vdata) == 1 )
                              vind = find( vdata );
                              method_params.voltage_timeseries = struct;
                              method_params.voltage_timeseries.value = h.data.tseries{vind}.name;
                              method_params.voltage_timeseries.name = 'voltage timeseries';
                              method_params.voltage_timeseries.type = 'string';
                           else
                              str = 'Set parameters to choose a voltage timeseries to match the AP templates to';
                              displayErrorMsg(str);
                              return;
                           end
                        end
                     else
                       method_params.voltage_timeseries = struct;
                       method_params.voltage_timeseries.value = h.data.tseries{end-1}.name;
                       method_params.voltage_timeseries.name  = 'voltage timeseries';
                       method_params.voltage_timeseries.type  = 'string';
                     end
                  end
                  voltage_name  = method_params.voltage_timeseries;
                  voltage_index = strcmpi(voltage_name.value, h.data.tseries_str);
                  if ~any(voltage_index)
                     str = 'Voltage timeseries not found, select a different timeseries';
                     displayErrorMsg(str);
                     return;
                  end
                  voltage_tseries     = h.data.tseries{voltage_index};
                  APtemplates         = tseries.data;
                  alignment_ind_templates = tseries.alignment_inds;
                  if isempty(APtemplates)
                     return;
                  end
                  str = sprintf( 'Extracting spikes using action potential templates\n' );
                  cprintf( printcol, str );
                  
                  normAPs  = tseries.params.normalise_aps.value; % get from AP generation
                  % get number of samples that each AP template is estimated from
                  APnumsamples = cellfun(@(f) size(f,2), tseries.APfamily);
                  
                  %copy template params over from APtemplate 
                  method_params.phasenumber = tseries.params.phasenumber;
                  method_params.first_phase_pos = tseries.params.first_phase_pos;
                  method_params.alignment_phase = tseries.params.alignment_phase;
                  
                  [APspikes, APstimes,alignment_ind_templates] = ...
                     extractSpikesUsingTemplates(APtemplates, APnumsamples,alignment_ind_templates, voltage_tseries, method, method_params, normAPs, h.options);
                  if isempty(APspikes)
                     return;
                  end
                  new_tseries.type    = 'spike';
                  new_tseries.data    = APspikes;
                  if iscell( APspikes{1} )
                     % if we ran outta memory we stored only spike times,
                     % which meant we lost the time vector, so I saved it in
                     % the fam struct
                     new_tseries.time = APspikes{1}{1}.time;
                  else
                     % gotta be cheeky getting time since we can plot spikes
                     % consecutively to make it easier to see how they change
                     for j = 1:size(APspikes,2)
                       if size(APspikes{j},1) > size(voltage_tseries.time,1)
                         APspikes{j} = APspikes{j}(1:size(voltage_tseries.time,1),1);
                       end
                     end
                     new_tseries.time = voltage_tseries.time( 1:size(APspikes{1}, 1) );
                  end
                  new_tseries.dt      = voltage_tseries.dt;
                  new_tseries.params  = method_params;
                  new_tseries.params.tool   = tool;
                  new_tseries.params.method = method;
                  new_tseries.APstimes= APstimes; % record spike times for getting rate later
                  tool_str            = [tseries.name '_spikes_' method];
                  tool_str            = [tseries.name '_spikes'];
                  instruct            = ['Creating ' tool_str ': rename?'];

               case 'merge templates'
                  str = sprintf( 'Merging action potential templates\n' );
                  cprintf( printcol, str );
                  % get time series to apply AP templates to
                  new_tseries = mergeAPtemplates( tseries, method, method_params );
                  if isempty(new_tseries)
                     return;
                  end
                  tool_str            = [tseries.name '_APs_' method];
                  tool_str            = [tseries.name '_APs'];
                  instruct            = ['Creating ' tool_str ': rename?'];

               case 'delete templates'
                  str = sprintf( 'Deleting action potential templates\n' );
                  cprintf( printcol, str );
                  % get time series to apply AP templates to
                  new_tseries = deleteAPtemplates( tseries, method, method_params );
                  if isempty(new_tseries)
                     return;
                  end
                  tool_str            = [tseries.name '_APs'];
                  instruct            = ['Creating ' tool_str ': rename?'];
               otherwise
                  displayErrorMsg('This tool don''t exist, give up');
                  return;
            end

         case 'spike'
            switch lower(tool)
               case 'firing rate'
                  str = sprintf( 'Calculating spiking rates\n' );
                  cprintf( printcol, str );
                  [rates, time, dt] = getRateFromSpikes(tseries, method, method_params);
                  if isempty(rates), return; end
                  new_tseries.type = 'rate';
                  new_tseries.data = rates;
                  new_tseries.time = time;
                  new_tseries.dt   = dt;
                  new_tseries.params= method_params;
                  new_tseries.params.tool   = tool;
                  new_tseries.params.method = method;
                  tool_str         = [tseries.name '_rates']; % _' method];
                  tool_str(tool_str=='_') = ' '; % make pretty for dialogue title
                  instruct         = ['Creating ' tool_str ': rename?'];

               case 'statistics'
                  str = sprintf( 'Calculating statistics\n' );
                  cprintf( printcol, str );
                  generateSpikeStatistics(tseries, method, method_params);
                  return;

               case 'spike operations'
                  str = sprintf( 'Applying spike operations\n' );
                  cprintf( printcol, str );
                  [spikes,stimes] = runSpikeOperations(tseries, method, method_params);

                  new_tseries.type    = 'spike';
                  new_tseries.data    = spikes;
                  new_tseries.time    = tseries.time;
                  new_tseries.dt      = tseries.dt;
                  new_tseries.params  = method_params;
                  new_tseries.params.tool   = tool;
                  new_tseries.params.method = method;
                  new_tseries.APstimes= stimes; % record spike times for getting rate later
                  tool_str            = [tseries.name '_spikes_' method];
                  tool_str            = [tseries.name '_spikes'];
                  instruct            = ['Creating ' tool_str ': rename?'];

              case 'export to excel'
                  str = sprintf( 'Exporting to excel\n' );
                  cprintf( printcol, str );
                  try
                     exportToExcel(tseries);
                  catch E
                     str = sprintf('\tAn unexpected error occurred while creating the excel file.\n');
                     caught_error = runtimeErrorHandler(E, 'message', str);
                     if ~isempty(caught_error), rethrow(caught_error); end
                  end
                  varargout = {h};
                  return;
               otherwise
            end

         case 'rate'
            str = sprintf( 'Generating rate statistics\n' );
            cprintf( printcol, str );
            % currently the only tool for 'rate' timeseries is generating
            % statistics, & the only statistic is autocorrelation
            generateRateStatistics(tseries, method, method_params);
            return;

         otherwise
            displayErrorMsg('You''re making shit up, this isn''t a data type');
            return;

      end
      elapsed = toc;fprintf('\tElapsed time: %0.4f seconds. Tool = %s\n',elapsed, tool);
      try
         if isfield(h.options, 'isBatch') && h.options.isBatch
           name = tool_str;
         else
           name = h.f.getFileName(instruct, tool_str, 60);
         end
      catch E
         if strcmp('MATLAB:inputdlg:InvalidInput',E.identifier)
            caught_error = runtimeErrorHandler(E,'ignore');
            if ~isempty(caught_error), rethrow(caught_error); end
            name = 'voltage';
         end
      end

      if isempty(name) % user cancelled process
         return
      end

      % if we're here we have a new timeseries - populate nec info & add to gui
      new_tseries.name = name;
      if isfield(tseries, 'params')
         new_tseries.params.nested_params = tseries.params;
      else
         new_tseries.params.nested_params = [];
      end

      %% save new tseries in user data & display in gui
      new_numtseries = h.data.num_tseries + 1;
      tseries_str    = h.data.tseries_str;
      used_names     = h.data.used_names;
      [name, used_names] = check4RepeatString(name, used_names);
      tseries_str{new_numtseries} = name;
      h.data.tseries{new_numtseries} = new_tseries;
      h.data.tseries_str = tseries_str;
      h.data.num_tseries = new_numtseries;
      h.data.used_names  = used_names;

      if strcmp('gui', h.f.uiType)
         set(h.curr_signal, 'String', tseries_str);
         set(h.curr_signal, 'Value',  new_numtseries);
         guidata(h.run_tool_button,h); % saves the change to handles
      else
         set(h.curr_signal, 'Items', tseries_str);
         set(h.curr_signal, 'Value',  tseries_str(new_numtseries));
      end
      h.f.curr_signal(h.curr_signal, [], h);

      if isfield(h.options, 'isBatch') && h.options.isBatch
        %guidata([], h);
        if strcmp('gui', h.f.uiType)
            guidata(h.run_tool_button, h); % saves the change to handles
        end
        varargout = {h};
      end
   end

   function save_voltage(h, varargin)
      % If gui:
      if nargin > 1
         hObject = varargin{1};
      end
      [tseries, ~, type, ts_name] = h.f.getCurrentVoltage(h);
      sname = title2Str(ts_name,1,1); % save name - options remove all punctuation
      if ~isfield(h.options, 'isBatch') || ~h.options.isBatch
        sname = h.f.getFileName('Name of saved variable in mat file ...', sname, 63);
      end

      if isempty(sname)
         return; % user's cancelled and hasn't provided a variable name
      end

      var_name = title2Str(ts_name,1,1,'_');
      % Check variable name size:
      if length(sname) > 50
         var_name = var_name(end - 50 : end);
         sname = sname(end - 50 : end);
      end
      eval_str = [sname ' = tseries;'];
      eval(eval_str);

      if ~isfield(h.options, 'isBatch') || ~h.options.isBatch
          displayErrorMsg( 'If you save into an existing smr file please ignore Matlab''s warning that it will be written over (select yes)' );

          if strcmpi(type,'voltage')
             filterspec = {'*.mat',  'MAT-files (*.mat)'; ...
                '*.smr',  'Spike files (*.smr)'; };
          else
             filterspec = {'*.mat'};
          end
          [fname, pname, findex] = uiputfile(filterspec, 'Save data as',...
             fullfile( h.data.last_dir, sname) );

          if isequal(fname,0) || isequal(pname,0) || findex==0 % (cancelled)
             return;
          end
        else
          % If its processing batches
          fname = sname;
          pname = h.options.batchPath;
          % Check for double backslashes
          dbl_bkslsh = strfind(pname,'\\');
          pname(dbl_bkslsh) = [];
          findex = 1;
        end

      full_name = fullfile(pname, fname);
      % Check if folder exists
      if ~exist(pname, 'dir')
         mkdir(pname);
      end

      % Check if file exists
      sufix = 0;
      valid_name = full_name;
      while exist([valid_name '.mat'], 'file')
         sufix = sufix + 1;
         valid_name = [full_name '_' num2str(sufix)];
      end
      full_name = valid_name;

      if findex==1 % .mat
         eval_str  = ['save(full_name, ''' sname ''', ''-v7.3'');'];
         eval(eval_str);

      elseif findex==2 % .smr
         % create a new empty file on disk
         dt = tseries.dt;
         scale = 1; offset = 0;

         % if file exists get first free channel to write to
         if exist( full_name, 'file' )
            [ fhand ]   = CEDS64Open( full_name, 0 ); % 0 for read/write mode (1 for read only)
            [ channel ] = CEDS64GetFreeChan( fhand );
            if fhand < 0
               displayErrorMsg( 'Error opening existing smr file %s (code %d)',...
                  fname, fhand );
               return;
            end

            % if file doesn't exist, create it & write to channel 1
         else
            channel = 1;
            [ fhand ] = CEDS64Create( full_name, 1, 0 );
            if fhand < 0
               displayErrorMsg( 'Error creating smr file %s (code %d)',...
                  fname, fillret );
               return;
            end
         end

         % set max time resolution of file
         dSecs         = CEDS64TimeBase( fhand, dt ); % sets file time base
         % set channel title & units
         [ok, units]   = CEDS64ChanTitle( fhand, channel, 'voltage' );
         % channel divider: 1/file_tick/channel_divider = sample_rate
         ok            = CEDS64SetWaveChan( fhand, channel, 1, 1, 1/dt );
         % get start time (kinda not nec coz we're not including an offset)
         starttime     = CEDS64SecsToTicks( fhand, 0 );
         % Convert to unit16: user_value = (channel_value) * scale /6553.6 + offset
         voltage       = int16( (tseries.data - offset) * 6553.6 / scale );
         % MUST write to ADC or real wave channel
         fillret       = CEDS64WriteWave( fhand, 1, voltage, starttime );
         ok            = CEDS64ChanComment( fhand, 1, 'Voltage generated in SpikeExtractionTool' );
         if fillret < 0
            CEDS64ErrorMessage(fillret);
            displayErrorMsg( 'Error writing to smr file %s (code %d)',...
               fname, fillret );
            return;
         end
         [ ok ] = CEDS64Close( fhand );
      end

      % Waveform channels hold data items that occur at fixed tick intervals,
      % - The sample interval is a fixed, integral multiple of the file tick.
      %   If the file tick were 1 microsecond, the available sample rates in
      %   the file would be 1000000/n where n is know as the channel divider.
      %   To achieve 500 Hz, with a 1 microsecond tick, the divider would be 2000.
      % - Currently, there are two types of waveform channels: Adc and RealWave.
      %   Adc channel type (code 1)
      %   These channels are designed to be efficient in data file space and store
      %   the data as 16-bit signed integers. You can set a scale and an offset
      %   value to convert these values into real, user units:
      %         user value = (16-bit value) * scale /6553.6 + offset
      %   With a scale of 1.0 and an offset of 0.0, the user values span the range
      %   -5.0 to 4.99985 user units.
      % - You can choose to read these channels as either 16-bit integer values or
      %   as 32-bit floating point values (converted with the scale and offset).
      % - You write these channels as 16-bit integers.
      % RealWave channel type (code 9)
      % - These channels store data as 32-bit floating point values.
      %   They also have a scale and offset that is used if you want to read the
      %   data back as 16-bit integers.
      %   16-bit value = (32-bit floating point value - offset) * 6553.6/scale
      % udpate last dir

      h.data.last_dir = pname;
      if strcmp('gui',h.f.uiType)
         guidata(hObject,h);
      end
   end

   function set_tool_params(h, varargin)
      if nargin > 1
         hObject = varargin{1};
      end

      % update curr tool now we're actually interested in this one
      if h.data.curr_tool ~= get(h.tool_list, 'Value')
         h.data.last_tool = h.data.curr_tool;
         h.data.curr_tool = get(h.tool_list, 'Value');
      end
      % for the particular choice of tool + implementation (i.e. method)
      % display the configuration parameters & allow user to set them
      [tseries, ts_num, data_type] = h.f.getCurrentVoltage(h);
      data_type    = tseries.type;

      if strcmp('gui', h.f.uiType)
         tool_num     = get(h.tool_list,'Value');
         tool_list    = get(h.tool_list,'String');
         tool         = tool_list{tool_num};

         method_num   = get(h.method_list,'Value');
         method_list  = get(h.method_list,'String');
         method       = method_list{method_num};
      else
         tool_list    = get(h.tool_list,'Items');
         tool         = get(h.tool_list,'Value');

         method_list  = get(h.method_list,'Items');
         method   = get(h.method_list,'Value');
      end

      params       = h.data.params; % current parameter values
   %   params = getDefaultToolParams; % allows to force to default params when testing
      method_params= getToolAndMethodParams(params, data_type, tool, method);
      dlg_name     = [tool ' using ' method];
        
      % if using AP templates for matched filtering need to select voltage
      % timseries to apply it to - add list to params for selection (unless
      % it's applied directly to a voltage timeseries)
      if strcmpi(tool, 'extract spikes')  && strcmpi(method, 'matched filter')
         if ~strcmpi(data_type, 'voltage')
            types     = getStructFieldFromCell(h.data.tseries, 'type');
            % names     = getStructFieldFromCell(h.data.tseries, 'name');
            % Previous line was causing an error when two voltages had same name
            % but different sufix.
            names = h.data.tseries_str;
            isvoltage = strcmpi('voltage', types);

            if sum(isvoltage)==0
               str    = 'No voltage time series to extract APs from, have another crack later';
               displayErrorMsg(str);
               return;
            end
            voltlist  = names(isvoltage); % list of all voltage time series'

            % create a parameter to ask user what voltage timeseries to apply AP templates to
            if ~isfield(method_params, 'voltage_timeseries')
               % find indices of all tseries that have type voltage
               voltage_timeseries.value    = names(find(isvoltage,1,'first'));
               voltage_timeseries.name     = 'voltage timeseries';
               voltage_timeseries.descript = 'select voltage timeseries to match APs to';
               voltage_timeseries.type     = 'list';
               voltage_timeseries.list     = names(isvoltage);
               voltage_timeseries.units    = 'timeseries name';
               method_params.voltage_timeseries = voltage_timeseries;

               % voltage timeseries parameter already exists - check it's still valid
            else
               voltage_timeseries = method_params.voltage_timeseries;
               % lists diff lengths so can't be equal (can't easily compare)
               if length(voltlist)~=length(voltage_timeseries.list)
                  voltage_timeseries.list = names(isvoltage);
                  voltage_tseries.value   = 1;
                  % lists same length but not same content - leave number choice the same
               elseif ~all(strcmpi(voltage_timeseries.list, voltlist))
                  voltage_timeseries.list = names(isvoltage);
               end
               method_params.voltage_timeseries = voltage_timeseries;
            end
         end
      elseif strcmpi(tool, 'export to excel')  && strcmpi(method, 'spike rate and count')
         displayErrorMsg('There are no parameters for this tool');
         return;
      end
      % Now we switch based on the tool and method.
      % some methods have special windows, others have simple inputdlg boxes
      switch tool
          % TODO: Fix this part of the code. Commented by Artemio 24 July 2024
          % case  'Identify AP templates'
          %     try 
          %         [method_params, cancel] = requestUserParamConfig_identify_AP_template(method_params, dlg_name);
          %     catch ME
          %       if strcmp('MATLAB:unassignedOutputs', ME.identifier)
          %         cancel = 1;
          %       else
          %         caught_error = runtimeErrorHandler(ME,'message', 'Something went wrong setting the parameters.');
          %         if ~isempty(caught_error), rethrow(caught_error); end
          %       end
          %     end    
          case 'extract_spikes'
              try 
                  
              catch ME
                if strcmp('MATLAB:unassignedOutputs', ME.identifier)
                  cancel = 1;
                else
                  caught_error = runtimeErrorHandler(ME,'message', 'Something went wrong setting the parameters.');
                  if ~isempty(caught_error), rethrow(caught_error); end
                end
              end    
          otherwise
              %these tools do not have special panels for them
              try
                [method_params, cancel] = requestUserParamConfig(method_params, dlg_name);
              catch ME
                if strcmp('MATLAB:unassignedOutputs', ME.identifier)
                  cancel = 1;
                else
                  cancel = 0; % Added by Artemio 24 July 2024  (trying to fix Identify IP templates param selection)
                  caught_error = runtimeErrorHandler(ME,'message', 'Something went wrong setting the parameters.');
                  if ~isempty(caught_error), rethrow(caught_error); end
                end
              end
      end
      
      

      % if merging templates, user must select template ID now we know how
      % many they want to merge (can't do simultaneously)
      names  = fieldnames( method_params );
      if ~cancel && any( strcmpi( 'number_of_templates_to_merge', names ) )
         mergeids = h.f.getUserTemplateMergeIDs( tseries, method_params );
         if isempty( mergeids )
            return;
         end
         % add merge ids to tseries datastruct since it isn't actually a user
         % parameter
         tseries.mergeAPs = mergeids;
         h.data.tseries{ ts_num } = tseries;
      elseif ~cancel && any( strcmpi( 'number_of_templates_to_remove', names ) )
         % If deleting templates
         deleteids = h.f.getUserTemplateDeleteIDs(tseries, method_params );
         if isempty( deleteids )
            return;
         end
         % add merge ids to tseries datastruct since it isn't actually a user
         % parameter
         tseries.deleteAPs = deleteids;
         h.data.tseries{ ts_num } = tseries;
      end

      params = setToolAndMethodParams(params, method_params, data_type, tool, method);
      h.data.params = params;
      if strcmp('gui', h.f.uiType)
         guidata(hObject,h); % saves changes to handles
      end
   end

   function time_min(h, varargin)
      if ~h.f.haveUserData(h)
         return;
      end

      if nargin > 1
         hObject = varargin{1};
      end

      % Get new max time & check with tseries time vector that it's within limits.
      prev_tlim = h.data.tlim;
      tseries   = h.f.getCurrentVoltage(h);
      tlim      = getTimeAndVoltageLimits(tseries, 'tlim');
      mint      = tlim(1); maxt = tlim(2);
      if strcmp('gui',h.f.uiType)
         str = get(h.time_min, 'String');
      else
         str = get(h.time_min, 'Value');
      end
      [ok, newt] = checkStringInput(str, 'float', mint, maxt);      
      if ~ok
         % users input is dodgy - gotta reverse engineer old text box min from
         % slider value & min display limit (zoom_min)
         % zoom_min = (text_max - text_min)*proport + text_min;
         % text_min = (zoom_min - text_max*proport) / (1-alpha)
         percent = get(h.time_slider, 'Value');
         proport = (1-percent)/2; % proportion max time changed by slider percentage
         zoom_min= h.data.tlim(1);
         if strcmp('gui', h.f.uiType)
            [~,text_max] = checkStringInput(get(h.time_max, 'String'), 'float');
            text_min= (zoom_min - proport*text_max) / (1-proport);
            displayErrorMsg(sprintf('Time must be between %d & %d', mint, maxt));
            set(h.time_min, 'String', sprintf('%.2f',text_min)); % reset to old val
         else
            [~,text_max] = checkStringInput(get(h.time_max, 'Value'), 'float');
            text_min= (zoom_min - proport*text_max) / (1-proport);
            displayErrorMsg(sprintf('Time must be between %d & %d', mint, maxt));
            set(h.time_min, 'Value', sprintf('%.2f',text_min)); % reset to old val
         end
         return
      end
      
      % if user's put in 0 assume they want the smallest time value, dt
      if strcmp(str,'0')
         newt = mint; ok = true;
         if strcmp('gui', h.f.uiType)
            set(h.time_min, 'String', sprintf('%.2f', tseries.dt));
         else
            set(h.time_min, 'Value', sprintf('%.2f', tseries.dt));
         end
      end
      
      % Gotta get min time from other text box because need to reset displayed
      % time lims just in case new min is larger than the last displayed max
      % (zoom allows you to zoom in from text box mins/maxes)
      if strcmp('gui', h.f.uiType)
         [~,maxt] = checkStringInput(get(h.time_max, 'String'), 'float'); % we know str is valid
         set(h.time_min, 'String', sprintf('%.2f',newt)); % reset to new val
         set(h.time_slider, 'Value', 1);
      else
         [~,maxt] = checkStringInput(get(h.time_max, 'Value'), 'float'); % we know str is valid
         set(h.time_min, 'Value', sprintf('%.2f',newt)); % reset to new val
         set(h.time_slider, 'Value', 1);
      end
      h.data.tlim(1) = newt;

      if ~(newt < maxt)
         currtseries = h.data.curr_tseries;
         maxt = h.data.tseries{currtseries}.time(end);
         if strcmp('gui', h.f.uiType)
           h.time_max.String = num2str(maxt);
         else
           set(h.time_max, 'Value', num2str(maxt));
         end
      end
      h.data.tlim(2) = maxt;

      % Update time sliders
      h = h.f.update_displacement_and_zoom(h, hObject);
      if strcmp('gui',h.f.uiType), guidata(hObject, h); end
      h = h.f.updateTimeSlider(h, hObject);
      h = updateSETFigure(h, tseries);

      % if there are any other SET gui's open, update their lims if current
      % time axes limits are the same for both (gui's own handle is always 1st)
      if strcmp('gui',h.f.uiType) && length(h.data.guihandles) > 1
         % remove any deleted handles
         invalid = ~ishandle( h.data.guihandles );
         h.data.guihandles( invalid ) = [];
         for ti=2:length(h.data.guihandles)
            other_gui     = h.data.guihandles(ti);
            other_handles = guidata(other_gui);
            other_data    = other_handles.data;
            other_tlim    = other_data.tlim;
            if all( compareFloats( prev_tlim, other_tlim, 0.1, 'percent' ) )
               set(other_h.time_min, 'String', sprintf('%.2f',newt)); % reset to new val
               set(other_h.time_slider, 'Value', 1);
               other_h.data.tlim(2) = maxt;
               other_h.data.tlim(1) = newt;
               othertseries    = h.f.getCurrentVoltage(other_handles);
               other_handles = updateSETFigure(other_handles, othertseries);
            end
         end
      end
   end

   function time_max(h,varargin)
   % "Maximum time" text box has changed, i.e. the user has changed the
   % zoom settings by manually entering a number in the time_max text box
      if nargin > 1
         hObject = varargin{1};
      end

      if ~h.f.haveUserData(h)
         return;
      end

      % Get new max time & check with tseries time vector that it's within limits.
      tseries    = h.f.getCurrentVoltage(h);
      prev_tlim  = h.data.tlim; % previous user time limits
      data_tlims = getTimeAndVoltageLimits(tseries, 'tlim'); % min/max poss time lims
      mint       = data_tlims(1); 
      maxt       = data_tlims(2);

      if strcmp('gui', h.f.uiType)
         str        = get(h.time_max, 'String');
      else
         str        = get(h.time_max, 'Value');
      end
      
      [ok, newt] = checkStringInput(str, 'float', mint, maxt); % Validate string
      
      if ~ok
         % users input is dodgy - gotta reverse engineer old text box max from
         % slider value & max display limit (zoom_max)
         % zoom_max   = text_max - (text_max - text_min)*proport;
         % zoom_max   = text_max - (text_max - text_min)*proport;
         percent = get(h.time_slider, 'Value');
         proport = (1-percent)/2; % proportion max time changed by slider percentage
         zoom_max= h.data.tlim(2);
         if strcmp('gui', h.f.uiType)
            [~,text_min] = checkStringInput(get(h.time_min, 'String'), 'float');
            text_max=  (zoom_max + proport*text_min) / (1-proport);% maxt;%text_max= (zoom_max + proport*text_min) / (1-proport);
            displayErrorMsg(sprintf('Time must be between %d & %d', mint, maxt));
            set(h.time_max, 'String', sprintf('%.2f',text_max)); % reset to old val
         else
            [~,text_min] = checkStringInput(get(h.time_min, 'Value'), 'float');
            text_max= maxt;% h.time_max.Value;%(zoom_max + proport*text_min) / (1-proport);
            displayErrorMsg(sprintf('Time must be between %d & %d', mint, maxt));
            set(h.time_max, 'Value', sprintf('%.2f',text_max)); % reset to old val
         end
         % Update time sliders
         h.f.updateTimeSlider(h);
         h = updateSETFigure(h, tseries);
         return
      end
      
      % Gotta get min time from other text box because need to reset displayed
      % time lims just in case new max is smaller than the last displayed min
      % (zoom allows you to zoom in from text box mins/maxes)
      if strcmp('gui', h.f.uiType)
         [~,mint] = checkStringInput(get(h.time_min, 'String'), 'float'); % we know str is valid
         set(h.time_max, 'String', sprintf('%.2f',newt)); % reset to new val
      else
         [~,mint] = checkStringInput(get(h.time_min, 'Value'), 'float'); % we know str is valid
         set(h.time_max, 'Value', sprintf('%.2f',newt)); % reset to new val
      end
      % set(h.time_slider, 'Value', 1);
      h.data.tlim(2) = newt;

      if ~(mint < newt)
         currtseries = h.data.curr_tseries;
         mint = h.data.tseries{currtseries}.time(1);
         if strcmp('gui', h.f.uiType)
           h.time_min.String = num2str(mint);
         else
           set(h.time_min, 'Value', num2str(mint));
         end
      end
      h.data.tlim(1) = mint;

      % Update time sliders
      h = h.f.update_displacement_and_zoom(h, hObject);
      if strcmp('gui', h.f.uiType), guidata(hObject, h); end
      h.f.updateTimeSlider(h);
      h = updateSETFigure(h, tseries);

      % if there are any other SET gui's open, update their lims if current
      % time axes limits are the same for both (gui's own handle is always 1st)
      if strcmp('gui', h.f.uiType) && length(h.data.guihandles) > 1
         invalid = ~ishandle( h.data.guihandles );
         h.data.guihandles( invalid ) = [];

         for ti=2:length( h.data.guihandles )
            other_gui     = h.data.guihandles(ti);
            % if other handle is valid update time if time lims match
            if ishandle(other_gui)
               other_handles = guidata(other_gui);
               other_data    = other_handles.data;
               other_tlim    = other_data.tlim;
               if all( compareFloats( prev_tlim, other_tlim, 0.1, 'percent' ) )
                  set(other_handles.time_max, 'String', sprintf('%.2f',newt)); % reset to new val
                  % time_max_Callback( other_handles.time_max, [], other_handles );

                  set(other_handles.time_max, 'String', sprintf('%.2f',newt)); % reset to new val
                  set(other_handles.time_slider, 'Value', 1);
                  other_handles.data.tlim(2) = newt;
                  other_handles.data.tlim(1) = mint;
                  othertseries    = h.f.getCurrentVoltage(other_handles);
                  other_handles = updateSETFigure(other_handles, othertseries);

               end

               % if we're here the handle has been deleted or something, so ditch
            else
               h.data.guihandles(ti) = [];
            end
         end
      end

      if strcmp('gui', h.f.uiType), guidata(hObject, h); end
   end

   function time_slider_updated(h, varargin)
      if nargin > 1
         hObject = varargin{1};
      end

      method = h.toggleZoomButton.UserData; % Loads current option (zoom or displacement)

      if strcmp('zoom',method)
         % slider zooms in/out to max/min values given in text boxes, so that
         % slider is a percentage of possible max/mins.
         percent = hObject.Value;
         if strcmp('gui',h.f.uiType)
            if percent == 1
               percent = 0.99999;
               h.time_slider.SliderStep(2) = 1;
               h.time_slider.SliderStep(1) = 5e-5;
               h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
            elseif percent > 0.999
               h.time_slider.SliderStep(2) = 1;
               h.time_slider.SliderStep(1) = 5e-5;
               h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
            elseif percent > 0.95
               h.time_slider.SliderStep(2) = 1;
               h.time_slider.SliderStep(1) = 1e-3;
               h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
            else
               h.time_slider.SliderStep = [0.01 0.1];
            end
         end
         min_zoom = 1e-5;
         h.data.zoomPercentage(1) = max(percent, min_zoom);%min(percent, min_zoom);

         prev_tlim  = h.data.tlim; % record time lims before slider was moved

         % calculate new min & max time lims by zoomin in from both ends
         tseries = h.f.getCurrentVoltage(h);
         if ~iscell(tseries.time)
            zoom_min = h.data.tlim(1);
            zoom_max = h.data.tlim(1) + ((tseries.time(end) - tseries.time(1))*(1-percent));

             h.data.tlim(1) = max( zoom_min, tseries.time(1) );
             h.data.tlim(2) = min( zoom_max, tseries.time(end) );
         end
         
      else
         % Displacement
         displacement = hObject.Value;
         prev_displacement = h.data.displacementPercentage(1);
         h.data.displacementPercentage(1) = displacement;
         prev_tlim  = h.data.tlim; % record time lims before slider was moved

         % calculate new min & max time lims by zoomin in from both ends
         tseries = h.f.getCurrentVoltage(h);
         disp_min = displacement * (tseries.time(end) - tseries.time(1));
         disp_max = disp_min + (h.data.tlim(2) - h.data.tlim(1));

         if disp_max > tseries.time(end)
            % Displacement has reached the maximum time.
            temp_max = disp_max;
            disp_max = tseries.time(end);
            disp_min = disp_max - (h.data.tlim(2) - h.data.tlim(1));
            
            if (prev_tlim(2) < disp_max) && h.data.time_slider_max
                % The slider has reached maximum, but user keeps pressing
                % button, so we will move the sliding bar to the end.
                h.data.time_slider_max = true;
                h.data.displacementPercentage(1) = 1; % Testing
                h.time_slider.Value = 1;
            elseif displacement < prev_displacement
                % The slider is at the max, but now the user is sliding to
                % the other direction (left). We need to make sure that the
                % slider actually moves.
                h.data.time_slider_max = false;
                disp_max = (displacement) * (tseries.time(end) - tseries.time(1));
                disp_min = disp_max - (h.data.tlim(2) - h.data.tlim(1));
            end
         end

         if disp_min < tseries.time(1)
            % Displacement has reached minimum time
            disp_min = tseries.time(1);
            disp_max = disp_min + (h.data.tlim(2) - h.data.tlim(1));
            
            if prev_tlim(1) ~= disp_min
                % The slider has reached maximum, but user keeps pressing
                % button, so we will move the sliding bar to the end.
                h.data.displacementPercentage(1) = 0; % Testing
                h.time_slider.Value = 0;
            end
         end

         h.data.tlim(1) = max( disp_min, tseries.time(1) );
         h.data.tlim(2) = min( disp_max, tseries.time(end) );

      end

      % Update time (horizontal) text boxes
      if strcmp('gui',h.f.uiType)
         h.time_min.String = sprintf('%0.2f',h.data.tlim(1));
         h.time_max.String = sprintf('%0.2f',h.data.tlim(2));
      else
         h.time_min.Value = sprintf('%0.2f',h.data.tlim(1));
         h.time_max.Value = sprintf('%0.2f',h.data.tlim(2));
      end

      tseries  = h.f.getCurrentVoltage(h);
      h        = updateSETFigure(h, tseries);
      this_gui = h.figure1;
      % if there are any other SET gui's open, update their lims if current
      % time axes limits are the same for both (gui's own handle is always 1st)
      if strcmp('gui', h.f.uiType) && length(h.data.guihandles) > 1
         for ti=2:length(h.data.guihandles)
            other_gui     = h.data.guihandles(ti);
            % if other handle is valid update time if time lims match
            if ishandle(other_gui) && ~isequal(this_gui, other_gui)
               other_handles = guidata(other_gui);
               other_data    = other_handles.data;
               other_tlim    = other_data.tlim;
               if all( compareFloats( prev_tlim, other_tlim, 0.1, 'percent' ) )
                  set( other_handles.time_slider, 'Value', percent );
                  try
                     other_handles.time_slider.Callback( other_handles.time_slider, other_handles );
%                      time_slider_Callback( other_handles.time_slider, [], other_handles );
                  catch ME
                     disp('check time slider update for other gui');
                  end
               end
               % if handle is invalid the gui's been deleted, so ditch
            else
               h.data.guihandles(ti) = [];
            end
         end
      end
      
      if strcmp('gui', h.f.uiType)
         guidata(hObject, h);
      end
   end

   function toggleZoomButton(h, hObject)
      % Toggles between arrows (displacement) and magnifier (zoom) icons
      warning('off','MATLAB:imagesci:png:libraryWarning'); % Ignore PNG associated warning

      % Get location of GUI files
      path = getFilePath();

      maxt = h.data.tlim(2);
      mint = h.data.tlim(1);
      zoom = 1 - h.data.zoomPercentage;

      if strcmp('zoom',hObject.UserData) % 
         % Displacement function has been selected
         if strcmp('gui', h.f.uiType)
            [x,~]=imread([path 'fig' filesep 'arrowsIcon.png']);% Load the displacement icon
            I2=imresize(x, [28.6 85.8]); % Resize icon
            hObject.CData = I2; % Assign icon to the button
            h.voltage_slider.SliderStep =  [0.01 0.1];
            
            % Check if the user has zoomed in, otherwise, there's nothing
            % to displace, i.e. the displacement step should be 0 and
            % the displacement bar should take the whole slider.
            if h.data.zoomPercentage(1) <= 1e-5
                % There's no zoom applied, i.e. the whole time series is
                % visible.
                h.time_slider.SliderStep = [0 1];
                h.data.displacementPercentage(1) = 0;
            else
                h.time_slider.SliderStep(2) = 1;
                h.time_slider.SliderStep(1) = zoom(1) / 3;% max(handles.time_slider.SliderStep(1) , handles.data.zoomPercentage(1));% Change the size of the vertical slider indicator to match the value zoomed in.
                h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*2);
            end
            % Update tooltip
            h = setTooltips(h, {hObject.Tag}, getTooltips({'toggleZoomButton_toZoom'}));
         else
            h.toggleZoomButton.Icon = [path 'fig' filesep 'arrowsIcon.png'];
            h = setTooltips(h, {'toggleZoomButton'}, getTooltips({'toggleZoomButton_toZoom'}));
         end
         hObject.UserData = 'disp'; % Change state to displacement
         h.time_slider.Value = h.data.displacementPercentage(1); % Update the position of the slider to represent displacement.
         h.voltage_slider.Value = 0.5; % Update the position of the slider to represent displacement.
         h.data.displacementPercentage(2) = 0;
         
         set(h.zoom_out_label,'Visible', 'off');
         set(h.zoom_in_label,'Visible', 'off');
      else
         if strcmp('gui', h.f.uiType)
            % Zoom function has been selected
            [x,~]=imread([path 'fig' filesep 'magnifierIcon.png']);% Load the displacement icon
            I2=imresize(x, [28.6 85.8]); % Resize icon
            hObject.CData = I2; % Assign icon to the button

            h.time_slider.SliderStep =  [0.001 0.1];%max(app.time_slider.SliderStep(1) , app.data.displacementPercentage(1)); % Change the size of the horizontal slider indicator to match the value displaced in.
            h.voltage_slider.SliderStep =  [0.01 0.1];%max(app.voltage_slider.SliderStep(1) , app.data.displacementPercentage(2)); %  Change the size of the vertical slider indicator to match the value displaced in.
         else
            h.toggleZoomButton.Icon = [path 'fig' filesep 'magnifierIcon.png'];
         end

         hObject.UserData = 'zoom'; % Change state to zoom
         h.time_slider.Value = h.data.zoomPercentage(1); % Update the position of the slider to represent zoom.
         h.voltage_slider.Value = h.data.zoomPercentage(2); % Update the position of the slider to represent zoom.
         
         set(h.zoom_out_label,'Visible', 'on');
         set(h.zoom_in_label,'Visible', 'on');
         
         if strcmp('gui', h.f.uiType)
           percent = 1 - zoom(1);
           if percent == 1
               percent = 0.99999;
               zoom(1) = 1 - percent;
               h.time_slider.SliderStep(2) = 1;
               h.time_slider.SliderStep(1) = 0.00005;
               h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
           elseif percent > 0.999
               h.time_slider.SliderStep(2) = 1;
               h.time_slider.SliderStep(1) = 0.00005;
               h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
           elseif percent > 0.95
               h.time_slider.SliderStep(2) = 1;
               h.time_slider.SliderStep(1) = 0.001;
               h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
           elseif percent > 0.9
               h.time_slider.SliderStep(2) = 1;
               h.time_slider.SliderStep(1) = 0.005;
               h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
           else
               h.time_slider.SliderStep = [0.01 0.1];
           end
           h.voltage_slider.SliderStep =  [0.01 0.1];
           
            % Update tooltip
            h = setTooltips(h, {hObject.Tag}, getTooltips({'toggleZoomButton_toDisplace'}));
         else
            h = setTooltips(h, {'toggleZoomButton'}, getTooltips({'toggleZoomButton_toDisplace'}));
         end
      end
   end

   function tool_list(hObject, h)
      % get name of tool chosen and update available methods
      [tseries, ~, data_type] = h.f.getCurrentVoltage(h);

      % Chec UI Type
      if strcmp('gui', h.f.uiType)
         tool_num  = get(hObject, 'Value');
         tool_name = get(hObject, 'String');
         % get methods list & update on gui, selecting first in list as currmethod
         methods   = getSETToolMethodsList(tool_name{tool_num}, data_type);
         set(h.method_list, 'String', methods);
         set(h.method_list, 'Value', 1);
      else
         tn = cellfun(@(x) strcmp(x,hObject.Value), hObject.Items);
         tool_num = find(tn == 1);
         tool_name = get(hObject, 'Items');
         % get methods list & update on gui, selecting first in list as currmethod
         methods   = getSETToolMethodsList(tool_name{tool_num}, data_type);
         set(h.method_list, 'Items', methods);
         set(h.method_list, 'Value', h.method_list.Items(1));
      end

      % update the tooltips
      tooltip_name = {['run_tool_' lower(tool_name{tool_num})]};
      h = setTooltips(h, {'run_tool_button'}, getTooltips(tooltip_name));
      tooltip_name = {['tool_list_' lower(tool_name{tool_num})]};
      h = setTooltips(h, {'tool_list'}, getTooltips(tooltip_name));

      guidata(hObject,h); % saves changes to handles
   end

   function h = updateTimeSlider(h, varargin)
      if nargin > 1
          hObject = varargin{1};
      end
      
      maxt    = h.data.tlim(2);
      mint    = h.data.tlim(1);
      trange  = maxt - mint;
      curr_ts = h.data.curr_tseries;

      h.data.zoomPercentage(1) = 1 - ((maxt-mint) / (h.data.tseries{curr_ts}.time(end) - h.data.tseries{curr_ts}.time(1)));
      h.data.displacementPercentage(1) = (maxt - trange) / h.data.tseries{curr_ts}.time(end);
      
      if strcmp( 'zoom', h.toggleZoomButton.UserData )
         h.time_slider.Value = h.data.zoomPercentage(1); % Update the position of the slider to represent zoom.

         percent = h.data.zoomPercentage(1);
         zoom = 1 - percent;
         if strcmp('gui', h.f.uiType) % Only if its GUI
           if percent == 1
             percent = 0.99999;
             zoom = 1 - percent;
             h.time_slider.SliderStep(2) = 1;
             h.time_slider.SliderStep(1) = 0.00005;
             h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
           elseif percent > 0.999
             h.time_slider.SliderStep(2) = 1;
             h.time_slider.SliderStep(1) = 0.00005;
             h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
           elseif percent > 0.9
             h.time_slider.SliderStep(2) = 1;
             h.time_slider.SliderStep(1) = 0.001;
             h.time_slider.SliderStep(2) = min(1, h.time_slider.SliderStep(1)*1.01);
           else
             h.time_slider.SliderStep = [0.01 0.1];
           end
         end
         
      else
        percent = h.data.zoomPercentage(1);
        zoom    = 1 - percent;
        
        if maxt >= h.data.tseries{curr_ts}.time(end-1)
          h.time_slider.Value = 1;
        elseif mint <= h.data.tseries{curr_ts}.time(1)
          h.time_slider.Value = 0;
        else
          h.time_slider.Value = h.data.displacementPercentage(1); % Update the position of the slider to represent displacement.
        end

        if strcmp( 'gui', h.f.uiType ) % Only if its GUI
          h.time_slider.SliderStep(2) = 1;
          h.time_slider.SliderStep(1) = zoom(1) / 3;% max(handles.time_slider.SliderStep(1) , handles.data.zoomPercentage(1));% Change the size of the vertical slider indicator to match the value zoomed in.
          h.time_slider.SliderStep(2) = min( 1, h.time_slider.SliderStep(1)*2 );
        end
      end
            
      if strcmp('gui', h.f.uiType) % Only if its GUI
        guidata(h.time_slider, h);
      end
      
      h.f.time_slider_updated(h, h.time_slider);
   end

   function updateVoltageSlider(h)
      maxv = h.data.vlim(2);
      minv = h.data.vlim(1);
      currtseries = h.data.curr_tseries;

      if strcmp( 'zoom', h.toggleZoomButton.UserData )
         try
            h.data.zoomPercentage(2) = 1 - ((maxv-minv) / (max(h.data.tseries{currtseries}.data) - min(h.data.tseries{currtseries}.data)));
            h.data.zoomPercentage(2) = max(h.data.zoomPercentage(2), 0); % Prevents the zoom to be less than 0
         catch E
            if strcmp('MATLAB:min:wrongInput', E.identifier) || strcmp('MATLAB:max:wrongInput', E.identifier)
               h.data.zoomPercentage(2) = 1 - ((maxv-minv) / (max(max(h.data.tseries{currtseries}.data{1}{1}.spikes)) - min(min(h.data.tseries{currtseries}.data{1}{1}.spikes))));
            else
               rethrow(E);
            end
         end
         h.voltage_slider.Value = h.data.zoomPercentage(2); % Update the position of the slider to represent zoom.
         if strcmp('gui', h.f.uiType), h.voltage_slider.SliderStep =  [0.1 0.2]; end%max(h.voltage_slider.SliderStep(1) , h.data.displacementPercentage(2)); %  Change the size of the vertical slider indicator to match the value displaced in.

      else
         try
            h.data.zoomPercentage(2) = 1 - ((maxv-minv) / (max(max(h.data.tseries{currtseries}.data{1}{1}.spikes)) - min(min(h.data.tseries{currtseries}.data{1}{1}.spikes))));
            h.data.displacementPercentage(2) = 1 - maxv / max(max(h.data.tseries{currtseries}.data{1}{1}.spikes));
         catch E
            if strcmp('MATLAB:max:wrongInput',E.identifier) || strcmp('MATLAB:min:wrongInput', E.identifier) || strcmp('MATLAB:cellRefFromNonCell',E.identifier)
               h.data.zoomPercentage(2) = 1 - ((maxv-minv) / (max(h.data.tseries{currtseries}.data) - min(h.data.tseries{currtseries}.data)));
               h.data.displacementPercentage(2) = 1 - maxv / max(h.data.tseries{currtseries}.data);
            else
               rethrow(E);
            end
         end
         h.voltage_slider.Value = h.data.displacementPercentage(2); % Update the position of the slider to represent displacement.
         if strcmp('gui', h.f.uiType)
            h.voltage_slider.SliderStep(2) = 0.2; % max(h.voltage_slider.SliderStep(1) , h.data.zoomPercentage(2));% Change the size of the vertical slider indicator to match the value zoomed in.
            h.voltage_slider.SliderStep(1) = 0.1; % * h.voltage_slider.SliderStep(2);
         end
      end
   end

   function voltage_slider_updated(h, hObject)
      method = h.toggleZoomButton.UserData; % Loads current option (zoom or displacement)
      
      doUpdateVoltageSlider = false;

      if strcmp('zoom',method)
         % slider zooms in/out to max/min values given in text boxes, so that
         % slider is a percentage of possible max/mins.
         percent = hObject.Value;
         h.data.zoomPercentage(2) = percent;
         if (0 < h.data.zoomPercentage(2)) && (0.1 > h.data.zoomPercentage(2))
             h.data.zoomPercentage(2) = 0.1;
             percent = 0.1;
             doUpdateVoltageSlider = true;
         end
         prev_vlim  = h.data.vlim; % record time lims before slider was moved

         % calculate new min & max time lims by zoomin in from both ends
         tseries = h.f.getCurrentVoltage(h);
         try
            vrange = max(tseries.data) - min(tseries.data);
            zoom_min = mean([h.data.vlim(1) h.data.vlim(2)]) - (vrange * (1-percent) / 2);
            zoom_max = mean([h.data.vlim(1) h.data.vlim(2)]) + (vrange * (1-min(percent,0.999)) / 2);

            h.data.vlim(1) = max( zoom_min, min(tseries.data));
            h.data.vlim(2) = min( zoom_max, max(tseries.data));
         catch E
            if strcmp('MATLAB:min:wrongInput', E.identifier) || strcmp('MATLAB:max:wrongInput', E.identifier)
               vrange = max(tseries.data{:},[],[1 2]) - min(tseries.data{:},[],[1 2]);
               zoom_min = mean([h.data.vlim(1) h.data.vlim(2)]) - (vrange * (1-percent) / 2);
               zoom_max = mean([h.data.vlim(1) h.data.vlim(2)]) + (vrange * (1-min(percent,0.999)) / 2);

               h.data.vlim(1) = max( zoom_min, min(tseries.data{:},[],[1 2]));
               h.data.vlim(2) = min( zoom_max, max(tseries.data{:},[],[1 2]));
            else
               rethrow(E);
            end
         end

         h.data.displacementPercentage(2) = 0;
      else
         % Displacement
         newdisplacement = hObject.Value;
         if newdisplacement > h.data.displacementPercentage(2)
            displacement = newdisplacement;%handles.data.displacementPercentage(2) - newdisplacement;
         else
            displacement = -newdisplacement;
         end
         h.data.displacementPercentage(2) = newdisplacement;
         prev_vlim  = h.data.vlim; % record time lims before slider was moved

         % calculate new min & max time lims by zoomin in from both ends
         tseries = h.f.getCurrentVoltage(h);
         disp_min = h.data.vlim(1) + displacement;
         disp_max = h.data.vlim(2) + displacement;
         if ~exist('vrange','var')
             vrange = disp_max - disp_min;
         end
            
         try
            h.data.vlim(1) = max( disp_min, min(tseries.data));
            h.data.vlim(2) = min( disp_max, max(tseries.data));
            if (h.data.vlim(2) < h.data.vlim(1))
                % If vlim(2) is smaller, it means vlim(1) was lower than
                % the limit, hence it was forced to be larger. 
                if h.data.vlim(1) == min(tseries.data)
                    % This line will leave vlim(1) as the bottom limit, and
                    % set vlim(2) to vlim(1) + vrange
                    h.data.vlim(2) = h.data.vlim(1) + vrange;
                else
                    % This line will leave vlim(2) as the top limit, and
                    % set vlim(1) to vlim(2) - vrange
                    h.data.vlim(2) = max(tseries.data);
                    h.data.vlim(1) = h.data.vlim(2) - vrange; 
                end
            end
         catch E
            if strcmp('MATLAB:min:wrongInput', E.identifier) || strcmp('MATLAB:max:wrongInput', E.identifier)
               h.data.vlim(1) = max( disp_min, min(tseries.data{:},[],[1 2]));
               h.data.vlim(2) = min( disp_max, max(tseries.data{:},[],[1 2]));
            else
               rethrow(E);
            end
         end
      end

      if strcmp('gui', h.f.uiType)
         % Update voltage (vertical) text boxes
         h.voltage_min.String = sprintf('%0.2f',h.data.vlim(1));
         h.voltage_max.String = sprintf('%0.2f',h.data.vlim(2));

         tseries = h.f.getCurrentVoltage(h);
         h = updateSETFigure(h, tseries);
         guidata(hObject, h);
      else
         % Update voltage (vertical) text boxes
         h.voltage_min.Value = sprintf('%0.2f',h.data.vlim(1));
         h.voltage_max.Value = sprintf('%0.2f',h.data.vlim(2));

         tseries = h.f.getCurrentVoltage(h);
         h = updateSETFigure(h, tseries);
      end

   end

   function voltage_max(h, hObject)
      if ~h.f.haveUserData(h)
         return;
      end

      % Get new max voltage & check with tseries data vector that it's within limits.
      tseries = h.f.getCurrentVoltage(h);
      if isnumeric( tseries.data )
         minv = min(tseries.data)*1.2;
         maxv = max(tseries.data)*1.2;
      elseif iscell( tseries.data )
         try
            minv = min( cellfun( @(d) min( d(:) ), tseries.data ) );
            maxv = max( cellfun( @(d) max( d(:) ), tseries.data ) );
         catch
            % If there is an error, it means the cell structure is a bit more
            % complicated. Set the min and max based on the first template
            minv = min(min(tseries.data{1}{1}.spikes)) * 1.2;
            maxv = max(max(tseries.data{1}{1}.spikes)) * 1.2;
         end
      end
      if strcmp('gui', h.f.uiType)
         str     = get(h.voltage_max, 'String');
      else
         str     = get(h.voltage_max, 'Value');
      end
      [ok, newv] = checkStringInput(str, 'float', minv, maxv);
      if ~ok
         % users input is dodgy - gotta reverse engineer old text box max from
         % slider value & max display limit (zoom_max)
         % zoom_max   = text_max - (text_max - text_min)*proport;
         % zoom_max   = text_max - (text_max - text_min)*proport;
         percent = get(h.voltage_slider, 'Value');
         proport = (1-percent)/2; % proportion max time changed by slider percentage
         zoom_max= h.data.vlim(2);
         if strcmp('gui', h.f.uiType)
            [~,text_min] = checkStringInput(get(h.voltage_min, 'String'), 'float');
            text_max = maxv; %(zoom_max + proport*text_min) / (1-proport);
            displayErrorMsg(sprintf('Voltage must be between %d & %d', minv, maxv));
            set(h.voltage_max, 'String', sprintf('%.2f',text_max)); % reset to old value
         else
            [~,text_min] = checkStringInput(get(h.voltage_min, 'Value'), 'float');
            text_max = maxv; %(zoom_max + proport*text_min) / (1-proport);
            displayErrorMsg(sprintf('Voltage must be between %d & %d', minv, maxv));
            set(h.voltage_max, 'Value', sprintf('%.2f',text_max)); % reset to old value
         end
         % Update voltage sliders
         h.f.updateVoltageSlider(h);
         h = updateSETFigure(h, tseries);
         if strcmp('gui', h.f.uiType), guidata(hObject, h); end
         return
      end
      % Gotta get min time from other text box because need to reset displayed
      % lims just in case new max is smaller than the last displayed min
      % (zoom allows you to zoom in from text box mins/maxes)
      if strcmp('gui', h.f.uiType)
         [~,minv] = checkStringInput(get(h.voltage_min, 'String'), 'float'); % we know str is valid
         set(h.voltage_max, 'String', sprintf('%.2f',newv)); % reset to new value
      end

      h.data.vlim(2) = newv;
      h.data.vlim(1) = minv;

      % Update voltage sliders
      h.f.updateVoltageSlider(h);
      
      h = updateSETFigure(h, tseries);
      if strcmp('gui', h.f.uiType), guidata(hObject, h); end
      
   end

   function voltage_min(h,hObject)
      if ~h.f.haveUserData(h)
         return;
      end

      % Get new max voltage & check with tseries data vector that it's within limits
      tseries = h.f.getCurrentVoltage(h);
      minv    = min(tseries.data)*1.2;
      maxv    = max(tseries.data)*1.2;

      if strcmp('gui', h.f.uiType)
         str     = get(h.voltage_min, 'String');
      else
         str     = get(h.voltage_min, 'Value');
      end

      [ok, newv] = checkStringInput(str, 'float', minv, maxv);
      % if user's put in 0 assume they want the smallest voltage value, dt
      if ~ok
         % users input is dodgy - gotta reverse engineer old text box min from
         % slider value & min display limit (zoom_min)
         % zoom_min = (text_max - text_min)*proport + text_min;
         % text_min = (zoom_min - text_max*proport) / (1-alpha)
         percent = get(h.voltage_slider, 'Value');
         proport = (1-percent)/2; % proportion max voltage changed by slider percentage
         zoom_min= h.data.vlim(1);
         if strcmp('gui', h.f.uiType)
            [~,text_max] = checkStringInput(get(h.voltage_max, 'String'), 'float');
            text_min = minv; % (zoom_min - proport*text_max) / (1-proport);
            displayErrorMsg(sprintf('Voltage must be between %d & %d', minv, maxv));
            set(h.voltage_min, 'String', sprintf('%.2f',text_min)); % reset to old val
         else
            [~,text_max] = checkStringInput(get(h.voltage_max, 'Value'), 'float');
            text_min = minv; % (zoom_min - proport*text_max) / (1-proport);
            displayErrorMsg(sprintf('Voltage must be between %d & %d', minv, maxv));
            set(h.voltage_min, 'Value', sprintf('%.2f',text_min)); % reset to old val
         end
         % Update voltage sliders
         h.f.updateVoltageSlider(h);
         h = updateSETFigure(h, tseries);
         if strcmp('gui', h.f.uiType), guidata(hObject, h); end
         return
      end
      % Gotta get max voltage from other text box because need to reset displayed
      % voltage lims just in case new max is smaller than the last displayed min
      % (zoom allows you to zoom in from text box mins/maxes)
      if strcmp('gui', h.f.uiType)
         [~,maxv] = checkStringInput(get(h.voltage_max, 'String'), 'float'); % we know str is valid
         set(h.voltage_min, 'String', sprintf('%.2f',newv)); % reset to new val
      else
         [~,maxv] = checkStringInput(get(h.voltage_max, 'Value'), 'float'); % we know str is valid
         set(h.voltage_min, 'Value', sprintf('%.2f',newv)); % reset to new val
      end

      h.data.vlim(1) = newv;
      if ~(newv < maxv)
         currtseries = h.data.curr_tseries;
         maxv = max(h.data.tseries{currtseries}.data);
         if strcmp('gui', h.f.uiType)
           h.voltage_max.String = sprintf('%.2f',maxv);
         else
           set(h.voltage_max, 'Value', sprintf('%.2f',maxv));
         end
      end
      h.data.vlim(2) = maxv;

      % Update voltage sliders
      h.f.updateVoltageSlider(h);

      h = updateSETFigure(h, tseries);
      if strcmp('gui', h.f.uiType), guidata(hObject, h); end
   end
   
   function h = update_displacement_and_zoom(h, hObject)
       % To update the values for horizontal displacement and zoom when the
       % user modifies the values in the textboxes
       maxt    = h.data.tlim(2);
       mint    = h.data.tlim(1);
       trange  = maxt - mint;
       curr_ts = h.data.curr_tseries;
       
       h.data.zoomPercentage(1) = 1 - ((maxt-mint) / (h.data.tseries{curr_ts}.time(end) - h.data.tseries{curr_ts}.time(1)));
       h.data.displacementPercentage(1) = (maxt - trange) / h.data.tseries{curr_ts}.time(end);
   end

end % Methods
end % ClassDef
