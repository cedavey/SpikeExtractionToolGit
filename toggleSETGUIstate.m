% Toggle the gui state from on to off and vice-versa. Off means it's in its
% initial state, before any data has been loaded. On means that we now have
% data and can display options.
% Inputs:
%  handles  - gui handles structure
%  state    - 'on' or 'off'
function handles = toggleSETGUIstate(handles,state)
%% UI objects
%  Select dataset
%    load_voltage    - press button to load new data file & erase current data
%    add_voltage     - button to add a new data file but retain all current data
%  Display voltage
%    curr_signal     - pop up list to choose current time series displayed
%    clear_voltage   - button to delete current voltage time series, but keep the rest
%  Generation parameters
%    param1(2,3,...) - text display of name of parameter used to generate display tseries
%    val1(2,3,...)   - text display of param value used to generate display tseries
%  Available tools
%    select_tool     - text display of select tool instruction
%    tool_list       - popup menu showing tools to apply to displayed tseries
%    select_method   - text display of select method instruction
%    method_list     - popup menu showing methods to apply current tool to displayed tseries
%    set_tool_params - button to launch gui to configure tool parameters
%  Figure
%    voltage_axes    - axes to plot time series in
%    voltage_slider  - slider to set axis limits for voltage
%    voltage_max     - text entry to allow user to set max voltage on axis
%    voltage_min     - text entry to allow user to set min voltage on axis
%    time_slider     - slider to set axis limits for time
%    time_max        - text entry to allow user to set max time on axis
%    time_min        - text entry to allow user to set min time on axis
%  AP templates
%    ap_title        - title for AP template family panel
%    curr_ap_family  - name of current AP template family to display

  % Check MATLAB's version
  v = ver('matlab');
  % Check if version is older than 2018b
  if str2double(v.Version) >= 9.5
      MLvers = 'current';
  else
      MLvers = 'old';
  end

   % Check if app or gui
   w = whos('handles');
   switch w.class
      case 'SpikeExtractionApp'
         uiType = 'app';
         % Check if data exists
         if handles.data.num_tseries>0
            data = handles.data;
            haveData = true;
         else
            haveData = false;
         end
      case 'SpikeExtractionTool_App'
         uiType = 'app';
         % Check if data exists
         if handles.data.num_tseries>0
            data = handles.data;
            haveData = true;
         else
            haveData = false;
         end
      case 'struct'
         uiType = 'gui';
         % Check if data exists
         if isfield(handles,'data') && handles.data.num_tseries>0
            data = handles.data;
            haveData = true;
         else
            haveData = false;
         end
      otherwise
         error('Couldn''t recognize the type of User Interface (gui or app)');
   end



   logic_state = strcmpi(state, 'on'); % logic is true if state is on, false otherwise

 % Select dataset
   set(handles.load_panel,     'Visible', 'on');
   set(handles.load_voltage_button,   'Visible', 'on'); % can always load a new datastruct
   set(handles.add_voltage_button,    'Visible', state);

 % Display voltage
   % curr_signal UI - setting depends on state we're toggling to
   set(handles.curr_signal,    'Visible', state);
   set(handles.clear_voltage_button,  'Visible', state);
   set(handles.voltage_panel,  'Visible', state);

 % Available tools
   set(handles.select_tool,    'Visible', state);
   set(handles.tool_list,      'Visible', state);
   set(handles.select_method,  'Visible', state);
   set(handles.method_list,    'Visible', state);
   set(handles.run_tool_button,       'Visible', state);
   set(handles.set_tool_params,'Visible', state);
   set(handles.tool_panel,     'Visible', state);

%  Figure
   set(handles.plot_panel,     'Visible', state);
   set(handles.voltage_slider, 'Visible', state);
   set(handles.voltage_max,    'Visible', state);
   set(handles.voltage_min,    'Visible', state);
   set(handles.time_slider,    'Visible', state);
   set(handles.time_max,       'Visible', state);
   set(handles.time_min,       'Visible', state);
   set(handles.save_voltage,   'Visible', state);
   set(handles.new_figure,     'Visible', state);
   set(handles.toggleZoomButton,'Visible', state);
   % if gui is deployed as stand alone then can't access variables
   if isdeployed
      set(handles.access_voltage, 'Visible', 'off');
   else
      set(handles.access_voltage, 'Visible', state);
   end
   set(handles.scroll_axes_slider,    'Visible', 'off');

   % The following panels don't get displayed until a voltage time series
   % is loaded or an AP template family is loaded - set to off until then
   % Generation parameters
   nP = 8;
   for i=1:8
      str = ['set(handles.param' num2str(1) ',''Visible'', ''off'');'];
      eval(str);
      str = ['set(handles.val' num2str(1)   ',''Visible'', ''off'');'];
      eval(str);
   end
   set(handles.param_panel,    'Visible', state);

   switch state
      case 'on'
         % Check Matlab version. Newer versions require true/false, whereas
         % older versions require 'on'/'off'
         if strcmp('current', MLvers), ctrl = true; else, ctrl = 'on'; end
         % Right click menu
         if strcmp('gui', uiType)
            handles.addVoltageMenu.Enable = ctrl;
            handles.clearVoltageMenu.Enable = ctrl;
         end
         % Menu Bar
         handles.addVoltageMenuBar.Enable = ctrl;
         handles.clearMenuBar.Enable = ctrl;
         handles.plotNewFigMenuBar.Enable = ctrl;
         handles.accessVoltageMenuBar.Enable = ctrl;
         handles.saveVoltageMenuBar.Enable = ctrl;
         % Slide bar update:
         handles.time_slider.Value = handles.data.zoomPercentage(1);
         handles.voltage_slider.Value = handles.data.zoomPercentage(2);

         % curr_signal UI
         if haveData
            % set font of drop down lists
            fontsize = 10;
            set(handles.tool_list,      'Fontsize', fontsize);
            set(handles.method_list,    'Fontsize', fontsize);
            set(handles.curr_signal,    'Fontsize', fontsize);

            % set font of text instruction boxes
            set(handles.select_tool,    'Fontsize', fontsize);
            set(handles.load_panel,     'Fontsize', fontsize);
            set(handles.load_voltage_button,   'Fontsize', fontsize);
            set(handles.add_voltage_button,    'Fontsize', fontsize);
            set(handles.clear_voltage_button,  'Fontsize', fontsize);
            set(handles.voltage_panel,  'Fontsize', fontsize);
            set(handles.select_tool,    'Fontsize', fontsize);
            set(handles.select_method,  'Fontsize', fontsize);
            set(handles.run_tool_button,       'Fontsize', fontsize);
            set(handles.set_tool_params,'Fontsize', fontsize);
            set(handles.save_voltage,   'Fontsize', fontsize);
            set(handles.new_figure,     'Fontsize', fontsize);
            set(handles.access_voltage, 'Fontsize', fontsize);

            % if have data, get type of current tseries & update tool lists
            data_type   = data.tseries{data.curr_tseries}.type;
            tool_list   = getSETToolList( data_type );
            method_list = getSETToolMethodsList(tool_list{1}, data_type);
            if strcmp('gui', uiType)
               set(handles.tool_list,  'String', tool_list);
               set(handles.tool_list,  'Value',  1);
               set(handles.method_list,'String', method_list);
               set(handles.method_list,'Value',  1);
               set(handles.voltage_min,'String', num2str(data.vlim(1)));
               set(handles.voltage_max,'String', num2str(data.vlim(2)));
               set(handles.time_min,   'String', num2str(data.tlim(1)));
               set(handles.time_max,   'String', num2str(data.tlim(2)));
            else
               set(handles.tool_list,  'Items', tool_list);
               set(handles.tool_list,  'Value',  tool_list(1));
               set(handles.method_list,'Items', method_list);
               set(handles.method_list,'Value',  method_list(1));
               set(handles.voltage_min,'Value', num2str(data.vlim(1)));
               set(handles.voltage_max,'Value', num2str(data.vlim(2)));
               set(handles.time_min,   'Value', num2str(data.tlim(1)));
               set(handles.time_max,   'Value', num2str(data.tlim(2)));
            end
         else
            if strcmp('gui', uiType)
               set(handles.tool_list,  'String',  'none');
               set(handles.tool_list,  'Value',   1);
               set(handles.method_list,'String',  'none');
               set(handles.method_list,'Value',   1);
            else
               set(handles.tool_list,  'Value',  handles.tool_list.Items(1));
               set(handles.method_list,  'Value',  handles.method_list.Items(1));
            end
         end

      case 'off'
         if strcmp('current', MLvers), ctrl = false; else, ctrl = 'off'; end
         
         if strcmp('gui', uiType)
           % Right click menu
           % Right click menu
           handles.addVoltageMenu.Enable = ctrl;
           handles.clearVoltageMenu.Enable = ctrl;
         end
         
         % Menu Bar
         handles.addVoltageMenuBar.Enable = ctrl;
         handles.clearMenuBar.Enable = ctrl;
         handles.plotNewFigMenuBar.Enable = ctrl;
         handles.accessVoltageMenuBar.Enable = ctrl;
         handles.saveVoltageMenuBar.Enable = ctrl;

         if strcmp('gui', uiType)
            set(handles.tool_list,  'String',  []);
            set(handles.tool_list,  'Value',   1);
            set(handles.method_list,'String',  []);
            set(handles.method_list,'Value',   1);
         else
            set(handles.tool_list,  'Items',  {''});
            set(handles.method_list,  'Items',  {''});
            set(handles.tool_list,  'Value',  handles.tool_list.Items(1));
            set(handles.method_list,  'Value',  handles.method_list.Items(1));
         end
   end

%  Generation parameters
%    param1(2,3,...) - text display of name of parameter used to generate display tseries
%    val1(2,3,...)   - text display of param value used to generate display tseries
%  Available tools
%    select_tool     - text display of select tool instruction
%    tool_list       - popup menu showing tools to apply to displayed tseries
%    select_method   - text display of select method instruction
%    method_list     - popup menu showing methods to apply current tool to displayed tseries
%    set_tool_params - button to launch gui to configure tool parameters
%  Figure
%    voltage_axes    - axes to plot time series in
%    voltage_slider  - slider to set axis limits for voltage
%    voltage_max     - text entry to allow user to set max voltage on axis
%    voltage_min     - text entry to allow user to set min voltage on axis
%    time_slider     - slider to set axis limits for time
%    time_max        - text entry to allow user to set max time on axis
%    time_min        - text entry to allow user to set min time on axis

    % if disabling GUI then get rid of all data in data
   if strcmpi(state,'off') && haveData
      %% Generate user data object to store
      data.num_tseries = 0;
      data.tseries     = [];   % initialise cell array of time series
      data.curr_tseries= [];   % current time series
      data.tseries_str = [];   % initialise names of time series
      data.used_names  = [];   % list of names already used by user
      data.last_tseries= [];   % last time series - restore if current ts deleted
      data.tlim        = [];   % time limits
      data.vlim        = [];   % voltage limits
      data.last_tool   = [];   % record of user's last tool
      data.curr_tool   = [];   % record of user's current tool
      data.ap_indices  = [];   % list of tseries that are actually ap templates
      data.params      = getDefaultToolParams; % default params for all tools & implementation methods

      if ~isfield(data,'last_dir')
         data.last_dir = pwd;  % keep if poss, else where they opened gui from
      end
      handles.data     = data; % record user data in handle

      % reset popupmenus
      if strcmp('gui', uiType)
         set(handles.curr_signal,'Visible','off');
         set(handles.curr_signal,'String', []);
         set(handles.curr_signal,'Value',  1);
         set(handles.tool_list,  'Visible','off');
         set(handles.tool_list,  'String', []);
         set(handles.tool_list,  'Value',  1);
         set(handles.method_list,'Visible','off');
         set(handles.method_list,'String', []);
         set(handles.method_list,'Value',  1);
      else

         set(handles.curr_signal,'Visible','off');
         handles.curr_signal.Items = {};
         handles.curr_signal.Value = {};
         set(handles.tool_list,  'Visible','off');
         handles.tool_list.Items = {};
         handles.tool_list.Value = {};
         set(handles.method_list,'Visible','off');
         handles.method_list.Items = {};
         handles.method_list.Value = {};
      end
   else % state=='on'
      if haveData
         tlim = handles.data.tlim;
         vlim = handles.data.vlim;

         % set time slider and time max/min text boxes from time length of data
         if strcmp('gui', uiType)
            set(handles.time_slider,    'Max',   1);
            set(handles.time_slider,    'Min',   0);
            set(handles.voltage_slider, 'Max',   1);
            set(handles.voltage_slider, 'Min',   0);
            set(handles.scroll_axes_slider,    'Max',   1);
            set(handles.scroll_axes_slider,    'Min',   0);
            set(handles.time_max,       'Value', tlim(2));
            set(handles.time_min,       'Value', tlim(1));
            set(handles.voltage_max,    'Value', vlim(2));
            set(handles.voltage_min,    'Value', vlim(1));
            set(handles.scroll_axes_slider,    'Value', 1);
         else
            handles.time_slider.Limits = [0 1];
            handles.voltage_slider.Limits = [0 1];
            handles.scroll_axes_slider.Limits = [0 1];
            set(handles.time_max,       'Value', num2str(tlim(2)));
            set(handles.time_min,       'Value', num2str(tlim(1)));
            set(handles.voltage_max,    'Value', num2str(vlim(2)));
            set(handles.voltage_min,    'Value', num2str(vlim(1)));
            set(handles.scroll_axes_slider,    'Value', 1);
         end

      end
   end

end
