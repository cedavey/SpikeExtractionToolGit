%% Questions
% - rescaling:
%       - some recordings have only +ve or -ve impacted by a rescaling
%         event. Does it make sense to rescale +ve & -ve separately?
%       - the system is designed to respond to slow changes, so can we
%         remove massive noise spikes prior, based on shape or crazy values
%         or something? Can we identify these events based on value or
%         might we get normalised voltage input or something?
% - positive before negative?
% - AP template:
%       - if you only plot the 1st 500 spikes contributing to the template,
%         & then plot the template on top, it looks really out --> the
%         templates are not stationary, but evolve over time
%
%% TO DO:
% - when extracting spikes, families can become v.similar in size so also
%   use correlation to determine which family within the template
% - also whene extracting spikes a one-off large spike spawns a new family,
%   which it perhaps shouldn't? Perhaps make allowed size diff. dependent
%   on noise?
% - button to reset axes (& perhaps to refresh?)
% - implement mean-shift clustering & PCA to separate neurons
% - if I find crazy voltage values, should I set them to 0?
% - allow user to remove AP templates from the set
% - align spikes based on min squared distance, rather than peak??
% - pointer in voltage plot window should give info on data
% - evaluate spike separation tool
% - enable parameter optimisation - will take ages to run!!
% - make a log file full of all the errors i've caught
% - when there's an error automatically email me the message & current status
% - test for empty results - APs, spikes
% - gotta change plotting spikes consecutively to a diff data type, rather
%   than 'spikes'
% - add button to link or unlink axes? Else only link with slider?
% - implement clustering using bayesian classifier, k-means, pca
% - should create objects so can pass pointer rather than timeseries
% - smooth firing rate estimates with gaussian kernel (the adaptive one!)
% - zoom & time range change simultaneously on multiple windows
% - check within AP family how correlated a spike is - choose largest!
% - enable displaying nested parameter record
% - enable putting red dots on spikes of timeseries for visualisation
% - allow resampling data to reduce dt & make things run faster
% - once a spike is detected it should be allocated to a template using a
%   priori probabilities as well as spike characteristic. The a priori
%   probabilities can be determined from the current spiking rate I guess?
%
%% Toolboxes required
% >> license('inuse')
% matlab
% statistics_toolbox
% wavelet_toolbox
% econometrics_toolbox
%
%% User data structure
%    data.num_tseries = 0;
%    data.tseries     = [];   % cell array of tseries with name, dt, method used to generate etc
%    data.curr_tseries= [];   % current time series as an index into tseries cell array
%    data.tseries_str = [];   % initialise names of time series
%    data.used_names  = [];   % list of names already used by user
%    data.last_tseries= [];   % last time series - restore if current ts deleted
%    data.tlim        = [];   % time limits
%    data.vlim        = [];   % voltage limits
%    data.last_tool   = [];   % record of user's last tool for voltage tseries (not aps)
%    data.curr_tool   = [];   % index of current tool
%    data.last_dir    = pwd;  % where they opened gui from
%    data.params      = [];   % record default or last chosen params for each tool/method combo
%    data.guihandles  = [];   % list of handles opened by gui
%    handles.data     = data; % record user data in handle
%
% Elements of the time series structure from which we get data for analysis
%    tseries.name     = name;
%    tseries.type     = 'voltage';
%    tseries.data     = voltage;
%    tseries.time     = time;
%    tseries.dt       = time(2) - time(1);
%    tseries.params   = [];


%% UI objects
%  Select dataset
%    load_voltage    - press button to load new data file & erase current data
%    add_voltage     - button to add a new data file but retain all current data
%  Display voltage
%    curr_signal     - pop up list to choose current time series displayed
%    clear_voltage   - button to delete current voltage time series, but keep the rest
%  Generation parameters
%    param1(2,3...)  - text display of name of parameter used to generate display tseries
%    val1(2,3,...)   - text display of param value used to generate display tseries
%  Available tools
%    select_tool     - text display of select tool instruction
%    tool_list       - popup menu showing tools to apply to displayed tseries
%    select_method   - text display of select method instruction
%    method_list     - popup menu showing methods to apply current tool to displayed tseries
%    run_tool        - initiate implementation of tool using select method
%    set_tool_params - button to launch gui to configure tool parameters
%  Figure
%    plot_panel      - panel within which we put a variable number of axes
%    voltage_slider  - slider to set axis limits for voltage
%    voltage_max     - text entry to allow user to set max voltage on axis
%    voltage_min     - text entry to allow user to set min voltage on axis
%    time_slider     - slider to set axis limits for time
%    time_max        - text entry to allow user to set max time on axis
%    time_min        - text entry to allow user to set min time on axis
%    save_voltage    - pops up a save dialogue box
%    scroll_axes     - if have more axes in figure than can be viewed allow scrolling
%    access_voltage  - saves current voltage to workspace
%    new_figure      - plot current tseries in new figure
%
% Notes:
% - to change callbacks for open, close etc, open gui in guide & select
% view -> View Callbacks, which shows all the figure's callback options
%
%

function varargout = SpikeExtractionTool(varargin)
% SPIKEEXTRACTIONTOOL MATLAB code for SpikeExtractionTool.fig
%      SPIKEEXTRACTIONTOOL, by itself, creates a new SPIKEEXTRACTIONTOOL or raises the existing
%      singleton*.
%
%      H = SPIKEEXTRACTIONTOOL returns the handle to a new SPIKEEXTRACTIONTOOL or the handle to
%      the existing singleton*.
%
%      SPIKEEXTRACTIONTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIKEEXTRACTIONTOOL.M with the given input arguments.
%
%      SPIKEEXTRACTIONTOOL('Property','Value',...) creates a new SPIKEEXTRACTIONTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpikeExtractionTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpikeExtractionTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpikeExtractionTool

% Last Modified by GUIDE v2.5 20-Jun-2019 11:42:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
   'gui_Singleton',  gui_Singleton, ...
   'gui_OpeningFcn', @SpikeExtractionTool_OpeningFcn, ...
   'gui_OutputFcn',  @SpikeExtractionTool_OutputFcn, ...
   'gui_LayoutFcn',  [] , ...
   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
   [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
   gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

% --- Executes just before SpikeExtractionTool is made visible.
function SpikeExtractionTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpikeExtractionTool (see VARARGIN)

% Choose default command line output for SpikeExtractionTool
handles.output   = hObject;

%% Generate user data object to store
data.num_tseries = 0;
data.tseries     = [];   % initialise cell array of time series
data.curr_tseries= [];   % current time series
data.tseries_str = [];   % initialise names of time series
data.used_names  = [];   % list of names already used by user
data.last_tseries= [];   % last time series - restore if current ts deleted
data.tlim        = [];   % time limits
data.vlim        = [];   % voltage limits
data.curr_tool   = [];   % record of user's current tool
data.last_tool   = [];   % record of user's last tool
data.params      = getDefaultToolParams; % default params for all tools & implementation methods
data.guihandles  = hObject;              % copy of all SET gui handles
data.zoomPercentage = [0 0]; % Records currently chosen zoom value
data.displacementPercentage = [0 0.5]; % Records currently chosen displacement value

data.last_dir    = pwd;  % where they opened gui from
handles.data     = data; % record user data in handle

handles = toggleSETGUIstate(handles,'off'); % switch everything off until data's loaded

% assign gui a name
set(hObject, 'Name', 'Spike Extraction Tool');
% set(gcf, 'menubar', 'figure') % turn the figure menu on (so can save etc) % I don't think it is a good idea to have this on, better to make a custom menu bar.


% prepare for user confirmation when GUI is closed
set(handles.figure1, 'CloseRequestFcn', @closeGUI);
%     set(gcf, 'WindowScrollWheelFcn', {@figure1_figScroll, handles});
removeToolBarButtons();

% from property inspector in guide, under WindowScrollWheelFcn
% @(hObject,eventdata)SpikeExtractionTool('figure1_WindowScrollWheelFcn',hObject,eventdata,guidata(hObject))

% Initialize options for debug and loading window
handles.loadingWindowOn = true;
handles.debugOption = 'semi';
handles.rescaleOption = 'at_end';

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpikeExtractionTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Starts logging everything to the log file: 'log_all.log'. 
% Everything means all the text printed in the command window, except for
% the one that is sent by the function 'printMessage'. Use such function
% when sending text to the command window that you don't want to be saved
% in the log file. i.e. instead of:
%                 cprintf('Keyword',str);
% try:
%                 printMessage('off','Keyword',str);
fid = fopen('./log_all.log', 'a'); % Opens log file to append this session's string
fprintf(fid, '\n\n-------------- %s @ %s | %s ---------------\n', getenv('Username'),getenv('UserDomain'),datestr(now, 0));
fclose(fid); % Close log file
diary('log_all.log'); % Activates the diary function, i.e. save all the activity into a file.
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
diary off % Deactivates the log. Hence nothing else will be added after closing the GUI.

% don't get user confirmation if there's no user data!
handles = guidata(src);
if ~isfield(handles, 'data_struct')
   delete(gcf);
   return;
end
if handles.data.num_tseries == 0
   delete(gcf);
   return;
end
response = userConfirmation('Do you want to close the GUI?',...
   'Close Request Function');
switch response
   case 'Yes'
      delete(gcf)
   case 'No'
      diary on % Reactivates the log
      return
end
end

function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
if ~haveUserData(handles), return; end

%     % see if pointer's in one of the axes
%     in_cp = isPointerInObject(handles, 'axes_coronal', 'crosshair');
%     if in_cp, return; end
%
%     % we're nowhere special so display pointer as an arrow
%     set(hObject,'pointer','arrow');
end

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
if ~haveUserData(handles), return; end

% see if pointer's in one of the axes
%     in_fig = isPointerInObject(handles, 'figure1');

%     in_lp = isPointerInObject(handles, 'axes_lateral');
%     if in_lp  % button press was in lateral axes
%         move_crosshairs(handles, 'lateral');
%         return
%     end

end

% Move crosshairs for each axes when user clicks left mouse button
% Because images are displayed with the starting point at the top LHS, I've
% flipped the data so that it starts from the bottom LHS. But now I need to
% keep the crosshairs flipped too, else voxel (1,1,1) is displayed as
% (maxX,maxY,maxZ), or something like that!
function [handles] = move_crosshairs(handles, axis)
% %     data_struct = handles.data_struct; % mem probs - avoid copying
tseries = getCurrentVoltage(handles);
maxt    = get(handles.time_slider,   'Max');
maxv    = get(handles.voltage_slider,'Max');

switch axis
   case 'lateral'
      %             % get new crosshair position values
      %             buttonpnt = get(handles.axes_lateral,'CurrentPoint');
      %             slice_Z   = round(buttonpnt(1,2)); if slice_Z==0, slice_Z=1; end
      %             slice_Y   = round(buttonpnt(1,1)); if slice_Y==0, slice_Y=1; end
      %
      %             % display new slice & update slice/slider values
      %             if (1 <= slice_Y && slice_Y <= max_Y) && ...
      %                     (1 <= slice_Z && slice_Z <= max_Z)
      %                 slice_X = round(str2num(get(handles.slice_X, 'String')));
      %                 % adjust for flipped images
      %                 slice_Z = max_Z+1-slice_Z;
      %                 str = sprintf('%3.3g', stat(slice_X, slice_Y, slice_Z));
      %                 set(handles.voxelValue, 'String', str);
      %                 updateImages(handles, stat, '', slice_Y, slice_Z);
      
   otherwise
      error('Unknown axes');
end % end switch axes
end


% checks if pointer is in object and, if so, returns true. If a pointer
% style is provided it also changes the pointer type.
% Input:
%   handles - figure handles
%   object  - name of object (i.e. a string)
%   pointer - a string containing new pointer type (optional)
function inobject = isPointerInObject(handles, object, pointer)
% E.g. position = get(handles.voxelValue_top, 'Position');
estr = ['pos = get(handles.' object ', ''Position'');'];
eval(estr);
pt   = get(handles.figure1, 'CurrentPoint'); % current pointer position

inobject = (pt(1,1)>=pos(1) && pt(1,1)<=(pos(1)+pos(3))) && ...
   (pt(1,2)>=pos(2) && pt(1,2)<=(pos(2)+pos(4)));
if inobject && exist('pointer', 'var')
   set(handles.figure1,'pointer',pointer);
end
end

function plot_panel_ButtonDownFcn(hObject, eventdata, handles)

end

% --- Outputs from this function are returned to the command line.
function varargout = SpikeExtractionTool_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.data;
varargout{2} = handles.output;
end

% --- Executes on button press in load_voltage.
function load_voltage_Callback(hObject, eventdata, handles)
% Warn if there is user data
okToClear = false;
if handles.data.num_tseries == 0
   okToClear = true;
end
% Ask for confirmation
if ~okToClear
   response = userConfirmation('All the currently loaded data will be cleared. Do you accept?', 'Don''t lose your data.');
   switch response
      case 'Yes'
         okToClear = true;
      case 'No'
         return
   end
end

mouseWaitingFunction(handles.figure1,@load_voltage,hObject,eventdata,handles);
end

function load_voltage(hObject, eventdata, handles)
last_dir = handles.data.last_dir; % keep before we write over it
old_numtseries = handles.data.num_tseries;

% get rid of all previous data
handles = toggleSETGUIstate(handles,'off');
handles.data.last_dir = last_dir;
guidata(hObject,handles);  % saves the change to handles
data    = handles.data;    % get user data from gui handle

% open smr, txt, or mat file
[data, success] = openVoltageFile(data);
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

guidata(hObject,handles);   % saves the change to handles
set(handles.curr_signal, 'String', data.tseries_str);
set(handles.curr_signal, 'Value',  1);
guidata(hObject,handles);

% curr_signal_Callback doesn't return handles so we have to save handles
% manually, then request a fresh copy using guidata
curr_signal_Callback(handles.curr_signal, '', handles);
guidata(hObject,handles);
end

% --- Executes on selection change in curr_signal.
function curr_signal_Callback(hObject, eventdata, handles)
if ~haveUserData(handles)
   return;
end
% find out what the user has selected to view
[tseries, ts_num, data_type] = getCurrentVoltage(handles);
% add current tseries to last viewed tseries before it's overwritten
handles.data.last_tseries    = handles.data.curr_tseries;
% update current tseries to tseries chosen by user
handles.data.curr_tseries    = ts_num;

% if data type of tseries hasn't changed then tool list doesn't need to
% change, but if we're displaying a new data type then the tool list and
% available methods for the tool needs to change
last_tseries = handles.data.tseries{handles.data.last_tseries};
last_type    = handles.data.tseries{handles.data.last_tseries}.type;
[last_tlim, last_vlim] = getTimeAndVoltageLimits(last_tseries);
[curr_tlim, curr_vlim] = getTimeAndVoltageLimits(tseries);
if ~strcmpi(data_type, last_type)
   % only remember last tool for voltage data cuz not many options for the others
   tool_num    = ternaryOp( strcmpi(data_type, 'voltage'), handles.data.last_tool, 1);
   if (isempty(tool_num) || tool_num==0), tool_num=1; end
   tool_list   = getSETToolList(data_type);
   tool        = tool_list{tool_num};
   method_list = getSETToolMethodsList(tool, data_type);
   set(handles.tool_list,   'String', tool_list);
   set(handles.tool_list,   'Value',  tool_num);
   set(handles.method_list, 'String', method_list);
   set(handles.method_list, 'Value', 1);
   
   %       % update min/max text boxes for voltage & time
   %       handles.data.tlim = curr_tlim;
   %       handles.data.vlim = curr_vlim;
   %       set(handles.voltage_min,'String', num2str(handles.data.vlim(1)));
   %       set(handles.voltage_max,'String', num2str(handles.data.vlim(2)));
   %       set(handles.time_min,   'String', num2str(handles.data.tlim(1)));
   %       set(handles.time_max,   'String', num2str(handles.data.tlim(2)));
   %
end

% update voltage & time sliders - keep them the same if timeseries length
% is the same as the last viewed tseries
if ~compareFloats(curr_tlim(2), last_tlim(2), 0.01, 'perc')
   % set time slider and time max/min text boxes from time length of data
   handles.data.tlim = curr_tlim;
   mint = handles.data.tlim(1); maxt = handles.data.tlim(2);
   set(handles.time_slider, 'Value',  0);
   set(handles.time_max,    'String', sprintf('%.2f',maxt));
   set(handles.time_min,    'String', sprintf('%.2f',mint));
end
if ~compareFloats(curr_vlim(2), last_vlim(2), 0.01, 'perc')
   % set time slider and time max/min text boxes from scale of data
   handles.data.vlim = curr_vlim;
   minv = handles.data.vlim(1); maxv = handles.data.vlim(2);
%    set(handles.voltage_slider, 'Value', 1);
   set(handles.voltage_max,    'String', sprintf('%.2f',maxv));
   set(handles.voltage_min,    'String', sprintf('%.2f',minv));
end
% udpate figure to show tseries chosen by user
guidata(hObject,handles); % saves changes to handles
handles = updateSETFigure(handles, tseries);
updateGUIParams(handles, tseries);

guidata(hObject,handles); % saves changes to handles
end

% --- Executes during object creation, after setting all properties.
function curr_signal_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', '');
set(hObject, 'Value', 1);
end

% --- Executes on button press in add_voltage.
function add_voltage_Callback(hObject, eventdata, handles)
mouseWaitingFunction(handles.figure1,@add_voltage,hObject,eventdata,handles);
end

function add_voltage(hObject,eventdata,handles)
last_tseries    = handles.data.curr_tseries;
new_num_tseries = handles.data.num_tseries + 1;

% open image files & retrieve matrices
[handles.data, success] = openVoltageFile(handles.data);

if success~=1
   str = sprintf('Error opening file, ignoring');
   displayErrorMsg(str);
   return
end
% loading new voltage tseries succeeded

% udpate drop down lists - gotta be done for success = 0 or 1 (?)
set(handles.curr_signal, 'String', handles.data.tseries_str);
% set voltage to 1st newly loaded voltage tseries
set(handles.curr_signal, 'Value', new_num_tseries);
curr_signal_Callback(handles.curr_signal, eventdata, handles);
end

% --- Executes on button press in clear_voltage.
function clear_voltage_Callback(hObject, eventdata, handles)
[tseries, ts_num, data_type, ts_name] = getCurrentVoltage(handles);
response = userConfirmation(['Delete ' ts_name '?'],...
   'Clear current time series?');
if strcmp(response,'No')
   return
end

mouseWaitingFunction(handles.figure1,@removeVoltage,handles); % Instead of removeVoltage(handles);

end

% --- Executes on selection change in tool_list.
function tool_list_Callback(hObject, eventdata, handles)
% get name of tool chosen and update available methods
[tseries, ~, data_type] = getCurrentVoltage(handles);
tool_num  = get(hObject, 'Value');
tool_name = get(hObject, 'String');
% get methods list & update on gui, selecting first in list as currmethod
methods   = getSETToolMethodsList(tool_name{tool_num}, data_type);
set(handles.method_list, 'String', methods);
set(handles.method_list, 'Value', 1);
guidata(hObject,handles); % saves changes to handles
end

% --- Executes during object creation, after setting all properties.
function tool_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in method_list.
function method_list_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function method_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in set_tool_params.
function set_tool_params_Callback(hObject, eventdata, handles)
% update curr tool now we're actually interested in this one
if handles.data.curr_tool ~= get(handles.tool_list, 'Value')
   handles.data.last_tool = handles.data.curr_tool;
   handles.data.curr_tool = get(handles.tool_list, 'Value');
end
% for the particular choice of tool + implementation (i.e. method)
% display the configuration parameters & allow user to set them
[tseries, ts_num, data_type] = getCurrentVoltage(handles);
data_type    = tseries.type;
tool_num     = get(handles.tool_list,'Value');
tool_list    = get(handles.tool_list,'String');
tool         = tool_list{tool_num};

method_num   = get(handles.method_list,'Value');
method_list  = get(handles.method_list,'String');
method       = method_list{method_num};

params       = handles.data.params; % current parameter values
% params = getDefaultToolParams;
method_params= getToolAndMethodParams(params, data_type, tool, method);
dlg_name     = [tool ' using ' method];

% if using AP templates for matched filtering need to select voltage
% timseries to apply it to - add list to params for selection (unless
% it's applied directly to a voltage timeseries)
if strcmpi(tool, 'extract spikes')  && strcmpi(method, 'matched filter')
   if ~strcmpi(data_type, 'voltage')
      types     = getStructFieldFromCell(handles.data.tseries, 'type');
      names     = getStructFieldFromCell(handles.data.tseries, 'name');
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
end

[method_params, cancel] = requestUserParamConfig(method_params, dlg_name);

% if merging templates, user must select template ID now we know how
% many they want to merge (can't do simultaneously)
names  = fieldnames( method_params );
if ~cancel && any( strcmpi( 'number_of_templates_to_merge', names ) )
   mergeids = getUserTemplateMergeIDs( tseries, method_params );
   if isempty( mergeids )
      return;
   end
   % add merge ids to tseries datastruct since it isn't actually a user
   % parameter
   tseries.mergeAPs = mergeids;
   handles.data.tseries{ ts_num } = tseries;
elseif ~cancel && any( strcmpi( 'number_of_templates_to_remove', names ) )
   % If deleting templates
   deleteids = getUserTemplateDeleteIDs( tseries, method_params );
   if isempty( deleteids )
      return;
   end
   % add merge ids to tseries datastruct since it isn't actually a user
   % parameter
   tseries.deleteAPs = deleteids;
   handles.data.tseries{ ts_num } = tseries;
end

params = setToolAndMethodParams(params, method_params, data_type, tool, method);
handles.data.params = params;
guidata(hObject,handles); % saves changes to handles
end

% if merging templates, user has chosen which template to merge with, &
% how many they'd like to merge - now ask the user which templates they
% want to actually merge with their chosen template (can't do this with
% the other params cuz need to know how many they want to merge before
% we can ask them which ones they want merged).
function userids = getUserTemplateMergeIDs( tseries, method_params )
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
   runtimeErrorHandler(ME,'ignore');
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

% if removing templates. User has chosen how many to remove.
function userids = getUserTemplateDeleteIDs( tseries, method_params )
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
   runtimeErrorHandler(ME,'ignore');
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

% --- Executes on button press in run_tool.
function run_tool_Callback(hObject, eventdata, handles)
mouseWaitingFunction(handles.figure1,@run_tool,hObject,eventdata,handles);
end

function run_tool(hObject, eventdata, handles)
% Implement the tool using the method & params requested
tool_list     = get(handles.tool_list, 'String');
tool_num      = get(handles.tool_list, 'Value');
tool          = tool_list{tool_num};

method_list   = get(handles.method_list, 'String');
method_num    = get(handles.method_list, 'Value');
method        = method_list{method_num};

[tseries, ~, type] = getCurrentVoltage(handles);
method_params = getToolAndMethodParams(handles.data.params, type, tool, method);

switch lower(type)
   case 'voltage'
      % only remember last tool for voltage coz others don't have many options
      handles.data.last_tool = handles.data.curr_tool;
      handles.data.curr_tool = get(handles.select_tool, 'Value');
      switch lower(tool)
         case 'rescale'
            [voltage, Rest] = rescaleVoltage(tseries, method, method_params, handles.debugOption, handles.loadingWindowOn);
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
            % get time series to apply AP templates to
            [APtemplates, APfamily] = identifyAPs(tseries, method, method_params,  handles.debugOption);
            if isempty(APtemplates)
               return;
            end
            new_tseries.type    = 'ap';
            new_tseries.data    = APtemplates;
            new_tseries.time    = (1:size(APtemplates,1))'*tseries.dt;
            new_tseries.dt      = tseries.dt;
            new_tseries.params  = method_params;
            new_tseries.params.tool   = tool;
            new_tseries.params.method = method;
            new_tseries.APfamily= APfamily;
            tool_str            = [tseries.name '_APs'];
            instruct            = ['Creating ' tool_str ': rename?'];
            
         case 'extract spikes'
            % allows user to generate spikes directly from voltage
            % timeseries, by generating AP templates enroute
            if strcmpi(method,'k means') 
               % Extracting spikes directly from voltage through K-means
               % clustering
               [APtemplates, APfamily]= identifyAPs(tseries, 'k means', method_params, handles.debugOption);
               if isempty(APtemplates)
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
               
            else
               [APtemplates, APfamily]= identifyAPs(tseries, 'threshold', method_params);
               if isempty(APtemplates)
                  return;
               end
               APnumsamples        = cellfun(@(f) size(f,2), APfamily);
               normAPs             = method_params.normalise_aps.value; % separate for when gen APs separately
               [APspikes,APtimes]  = extractSpikesUsingTemplates( APtemplates, APnumsamples, tseries, method, method_params, normAPs );
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
               str = 'Set parameters to choose a voltage timeseries to match the AP templates to';
               displayErrorMsg(str);
               return;
            end
            voltage_name        = method_params.voltage_timeseries;
            voltage_index       = strcmpi(voltage_name.value, handles.data.tseries_str);
            if ~any(voltage_index)
               str = 'Voltage timeseries not found, select a different timeseries';
               displayErrorMsg(str);
               return;
            end
            voltage_tseries     = handles.data.tseries{voltage_index};
            APtemplates         = tseries.data;
            if isempty(APtemplates)
               return;
            end
            normAPs             = tseries.params.normalise_aps.value; % get from AP generation
            % get number of samples that each AP template is estimated from
            APnumsamples        = cellfun(@(f) size(f,2), tseries.APfamily);
            [APspikes, APstimes]= ...
               extractSpikesUsingTemplates(APtemplates, APnumsamples, voltage_tseries, method, method_params, normAPs);
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
            % get time series to apply AP templates to
            new_tseries = mergeAPtemplates( tseries, method, method_params );
            if isempty(new_tseries)
               return;
            end
            tool_str            = [tseries.name '_APs_' method];
            tool_str            = [tseries.name '_APs'];
            instruct            = ['Creating ' tool_str ': rename?'];
         
         case 'delete templates'
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
            generateSpikeStatistics(tseries, method, method_params);
            return;
            
         case 'spike operations'
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
            
         otherwise
      end
      
   case 'rate'
      % currently the only tool for 'rate' timeseries is generating
      % statistics, & the only statistic is autocorrelation
      generateRateStatistics(tseries, method, method_params);
      return;
      
   otherwise
      displayErrorMsg('You''re making shit up, this isn''t a data type');
      return;
      
end

name = getFileName(instruct, tool_str, 60);
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
new_numtseries = handles.data.num_tseries + 1;
tseries_str    = handles.data.tseries_str;
used_names     = handles.data.used_names;
[name, used_names] = check4RepeatString(name, used_names);
tseries_str{new_numtseries} = name;
handles.data.tseries{new_numtseries} = new_tseries;
handles.data.tseries_str = tseries_str;
handles.data.num_tseries = new_numtseries;
handles.data.used_names  = used_names;

set(handles.curr_signal, 'String', tseries_str);
set(handles.curr_signal, 'Value',  new_numtseries);
guidata(handles.run_tool,handles); % saves the change to handles
curr_signal_Callback(handles.curr_signal, [], handles);

end

% --- Executes on button press in save_voltage.
function save_voltage_Callback(hObject, eventdata, handles)
[tseries, ~, type, ts_name] = getCurrentVoltage(handles);
sname = title2Str(ts_name,1,1); % save name - options remove all punctuation
sname = getFileName('Name of saved variable in mat file ...', sname, 63);
if isempty(sname)
   return; % user's cancelled and hasn't provided a variable name
end

var_name = title2Str(ts_name,1,1,'_');
eval_str = [sname ' = tseries;'];
eval(eval_str);

displayErrorMsg( 'If you save into an existing smr file please ignore Matlab''s warning that it will be written over (select yes)' );

if strcmpi(type,'voltage')
   filterspec = {'*.mat',  'MAT-files (*.mat)'; ...
      '*.smr',  'Spike files (*.smr)'; };
else
   filterspec = {'*.mat'};
end
[fname, pname, findex] = uiputfile(filterspec, 'Save data as',...
   fullfile( handles.data.last_dir, sname) );

if isequal(fname,0) || isequal(pname,0) || findex==0 % (cancelled)
   return;
end
full_name = fullfile(pname, fname);
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

handles.data.last_dir = pname;
guidata(hObject,handles);

end

% --- Executes on slider movement.
function time_slider_Callback(hObject, eventdata, handles)
if ~haveUserData(handles)
   return;
end
mouseWaitingFunction(handles.figure1,@time_slider_updated,hObject,eventdata,handles);
end

function time_slider_updated(hObject,eventdata,handles)
method = handles.toggleZoomButton.UserData; % Loads current option (zoom or displacement)

if strcmp('zoom',method)
   % slider zooms in/out to max/min values given in text boxes, so that
   % slider is a percentage of possible max/mins.
   percent = hObject.Value;
   handles.data.zoomPercentage(1) = percent;
   prev_tlim  = handles.data.tlim; % record time lims before slider was moved
   
   % calculate new min & max time lims by zoomin in from both ends   
   tseries = getCurrentVoltage(handles);
   zoom_min = handles.data.tlim(1);
   zoom_max = handles.data.tlim(1) + ((tseries.time(end) - tseries.time(1))*(1-percent));
   
%    if zoom_max > tseries.time(end)
%       zoom_min = tseries.time(end) - zoom_max;
%       zoom_max = tseries.time(end);
%    end
%    
%    if zoom_min < tseries.time(1)
%       zoom_max = tseries.time(1) + zoom_min;
%       zoom_min = tseries.time(1);
%    end
   
   handles.data.tlim(1) = max( zoom_min, tseries.time(1) );
   handles.data.tlim(2) = min( zoom_max, tseries.time(end) );
else
   % Displacement
   displacement = hObject.Value;
   handles.data.displacementPercentage(1) = displacement;
   prev_tlim  = handles.data.tlim; % record time lims before slider was moved
   
   % calculate new min & max time lims by zoomin in from both ends   
   tseries = getCurrentVoltage(handles);
   disp_min = displacement * (tseries.time(end) - tseries.time(1));
   disp_max = disp_min + (handles.data.tlim(2) - handles.data.tlim(1));
   
   if disp_max > tseries.time(end)
      disp_max = tseries.time(end);
      disp_min = disp_max - (handles.data.tlim(2) - handles.data.tlim(1));
   end
   
   if disp_min < tseries.time(1)
      disp_min = tseries.time(1);
      disp_max = disp_min + (handles.data.tlim(2) - handles.data.tlim(1));
   end
   
   handles.data.tlim(1) = max( disp_min, tseries.time(1) );
   handles.data.tlim(2) = min( disp_max, tseries.time(end) );
     
end
tseries = getCurrentVoltage(handles);
handles = updateSETFigure(handles, tseries);
   
% if there are any other SET gui's open, update their lims if current
% time axes limits are the same for both (gui's own handle is always 1st)
if length(handles.data.guihandles) > 1
   for ti=2:length(handles.data.guihandles)
      other_gui     = handles.data.guihandles(ti);
      % if other handle is valid update time if time lims match
      if ishandle(other_gui)
         other_handles = guidata(other_gui);
         other_data    = other_handles.data;
         other_tlim    = other_data.tlim;
         if all( compareFloats( prev_tlim, other_tlim, 0.1, 'percent' ) )
            set( other_handles.time_slider, 'Value', percent );
            time_slider_Callback( other_handles.time_slider, [], other_handles );
         end
         % if handle is invalid the gui's been deleted, so ditch
      else
         handles.data.guihandles(ti) = [];
      end
   end
end

guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function time_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function time_max_Callback(hObject, eventdata, handles)
if ~haveUserData(handles)
   return;
end

% Get new max time & check with tseries time vector that it's within limits.
tseries    = getCurrentVoltage(handles);
prev_tlim  = handles.data.tlim; % previous user time limits
data_tlims = getTimeAndVoltageLimits(tseries, 'tlim'); % min/max poss time lims
mint       = data_tlims(1); maxt = data_tlims(2);
str        = get(handles.time_max, 'String');
[ok, newt] = checkStringInput(str, 'float', mint, maxt);
if ~ok
   % users input is dodgy - gotta reverse engineer old text box max from
   % slider value & max display limit (zoom_max)
   % zoom_max   = text_max - (text_max - text_min)*proport;
   % zoom_max   = text_max - (text_max - text_min)*proport;
   percent = get(handles.time_slider, 'Value');
   proport = (1-percent)/2; % proportion max time changed by slider percentage
   zoom_max= handles.data.tlim(2);
   [~,text_min] = checkStringInput(get(handles.time_min, 'String'), 'float');
   text_max= (zoom_max + proport*text_min) / (1-proport);
   displayErrorMsg(sprintf('Time must be between %d & %d', mint, maxt));
   set(handles.time_max, 'String', sprintf('%.2f',text_max)); % reset to old val
   return
end
% Gotta get min time from other text box because need to reset displayed
% time lims just in case new max is smaller than the last displayed min
% (zoom allows you to zoom in from text box mins/maxes)
[~,mint] = checkStringInput(get(handles.time_min, 'String'), 'float'); % we know str is valid
set(handles.time_max, 'String', sprintf('%.2f',newt)); % reset to new val
% set(handles.time_slider, 'Value', 1);
handles.data.tlim(2) = newt;
handles.data.tlim(1) = mint;
handles = updateSETFigure(handles, tseries);

% if there are any other SET gui's open, update their lims if current
% time axes limits are the same for both (gui's own handle is always 1st)
if length(handles.data.guihandles) > 1
   invalid = ~ishandle( handles.data.guihandles );
   handles.data.guihandles( invalid ) = [];
   
   for ti=2:length( handles.data.guihandles )
      other_gui     = handles.data.guihandles(ti);
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
            othertseries    = getCurrentVoltage(other_handles);
            other_handles = updateSETFigure(other_handles, othertseries);
            
         end
         
         % if we're here the handle has been deleted or something, so ditch
      else
         handles.data.guihandles(ti) = [];
      end
   end
end

guidata(hObject, handles);
end

function time_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
end

function time_min_Callback(hObject, eventdata, handles)
if ~haveUserData(handles)
   return;
end

% Get new max time & check with tseries time vector that it's within limits.
prev_tlim = handles.data.tlim;
tseries   = getCurrentVoltage(handles);
tlim      = getTimeAndVoltageLimits(tseries, 'tlim');
mint      = tlim(1); maxt = tlim(2);
str       = get(handles.time_min, 'String');
[ok, newt] = checkStringInput(str, 'float', mint, maxt);
% if user's put in 0 assume they want the smallest time value, dt
if strcmp(str,'0')
   newt = mint; ok = true;
   set(handles.time_min, 'String', sprintf('%.2f', tseries.dt));
end
if ~ok
   % users input is dodgy - gotta reverse engineer old text box min from
   % slider value & min display limit (zoom_min)
   % zoom_min = (text_max - text_min)*proport + text_min;
   % text_min = (zoom_min - text_max*proport) / (1-alpha)
   percent = get(handles.time_slider, 'Value');
   proport = (1-percent)/2; % proportion max time changed by slider percentage
   zoom_min= handles.data.tlim(1);
   [~,text_max] = checkStringInput(get(handles.time_max, 'String'), 'float');
   text_min= (zoom_min - proport*text_max) / (1-proport);
   displayErrorMsg(sprintf('Time must be between %d & %d', mint, maxt));
   set(handles.time_min, 'String', sprintf('%.2f',text_min)); % reset to old val
   return
end
% Gotta get min time from other text box because need to reset displayed
% time lims just in case new min is larger than the last displayed max
% (zoom allows you to zoom in from text box mins/maxes)
[~,maxt] = checkStringInput(get(handles.time_max, 'String'), 'float'); % we know str is valid
set(handles.time_min, 'String', sprintf('%.2f',newt)); % reset to new val
set(handles.time_slider, 'Value', 1);
handles.data.tlim(1) = newt;
handles.data.tlim(2) = maxt;
handles = updateSETFigure(handles, tseries);
guidata(hObject, handles);

% if there are any other SET gui's open, update their lims if current
% time axes limits are the same for both (gui's own handle is always 1st)
if length(handles.data.guihandles) > 1
   % remove any deleted handles
   invalid = ~ishandle( handles.data.guihandles );
   handles.data.guihandles( invalid ) = [];
   for ti=2:length(handles.data.guihandles)
      other_gui     = handles.data.guihandles(ti);
      other_handles = guidata(other_gui);
      other_data    = other_handles.data;
      other_tlim    = other_data.tlim;
      if all( compareFloats( prev_tlim, other_tlim, 0.1, 'percent' ) )
         set(other_handles.time_min, 'String', sprintf('%.2f',newt)); % reset to new val
         set(other_handles.time_slider, 'Value', 1);
         other_handles.data.tlim(2) = maxt;
         other_handles.data.tlim(1) = newt;
         othertseries    = getCurrentVoltage(other_handles);
         other_handles = updateSETFigure(other_handles, othertseries);
      end
   end
end
end


% --- Executes during object creation, after setting all properties.
function time_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
end

% --- Executes on slider movement.
function voltage_slider_Callback(hObject, eventdata, handles)
if ~haveUserData(handles)
   return;
end

mouseWaitingFunction(handles.figure1,@voltage_slider_updated,hObject,eventdata,handles);
% % slider zooms in/out to max/min values given in text boxes, so that
% % slider is a percentage of possible max/mins. If move slider to 80%,
% % then increase min by 10% and reduce max by 10%
% percent    = get(handles.voltage_slider, 'Value'); % get old voltage
% str        = get(handles.voltage_min, 'String');
% [ok, minv] = checkStringInput(str, 'float');
% str        = get(handles.voltage_max, 'String');
% [ok, maxv] = checkStringInput(str, 'float');
% 
% tseries    = getCurrentVoltage( handles );
% switch lower( tseries.type )
%    case 'voltage'
%       zoom_min   =  minv * percent;
%       zoom_max   =  maxv * percent;
%    otherwise
%       zoom_min   = (maxv - minv)*((1-percent)/2) + minv;
%       zoom_max   =  maxv - (maxv - minv)*((1-percent)/2);
% end
% if iscell( tseries.data )
%    % get min/max from 1st column of each cell's matrix, cuz it's taking
%    % f*cking ages t
%    handles.data.vlim(1) = max( zoom_min, min( cellfun(@(d) min(d(:,1)), tseries.data ) ) );
%    handles.data.vlim(2) = min( zoom_max, max( cellfun(@(d) max(d(:,1)), tseries.data ) ) );
% else
%    handles.data.vlim(1) = max( zoom_min, min(tseries.data(:)) );
%    handles.data.vlim(2) = min( zoom_max, max(tseries.data(:)) );
% end
end

function voltage_slider_updated(hObject,eventdata,handles)
method = handles.toggleZoomButton.UserData; % Loads current option (zoom or displacement)

if strcmp('zoom',method)
   % slider zooms in/out to max/min values given in text boxes, so that
   % slider is a percentage of possible max/mins.
   percent = hObject.Value;
   handles.data.zoomPercentage(2) = percent;
   prev_vlim  = handles.data.vlim; % record time lims before slider was moved
   
   % calculate new min & max time lims by zoomin in from both ends   
   tseries = getCurrentVoltage(handles);
   vrange = max(tseries.data) - min(tseries.data);
   zoom_min = mean([handles.data.vlim(1) handles.data.vlim(2)]) - vrange * (1-percent);
   zoom_max = mean([handles.data.vlim(1) handles.data.vlim(2)]) + vrange * (1-percent);

   handles.data.vlim(1) = max( zoom_min, min(tseries.data));
   handles.data.vlim(2) = min( zoom_max, max(tseries.data));
else
   % Displacement
   newdisplacement = hObject.Value;
   displacement = handles.data.displacementPercentage(2) - newdisplacement;
   handles.data.displacementPercentage(2) = newdisplacement;
   prev_vlim  = handles.data.vlim; % record time lims before slider was moved
   
   % calculate new min & max time lims by zoomin in from both ends   
   tseries = getCurrentVoltage(handles);
   disp_min = handles.data.vlim(1) + displacement;
   disp_max = handles.data.vlim(2) + displacement;
      
   handles.data.vlim(1) = max( disp_min, min(tseries.data));
   handles.data.vlim(2) = min( disp_max, max(tseries.data));
end

tseries = getCurrentVoltage(handles);
handles = updateSETFigure(handles, tseries);
guidata(hObject, handles);

end

% --- Executes during object creation, after setting all properties.
function voltage_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function voltage_max_Callback(hObject, eventdata, handles)
if ~haveUserData(handles)
   return;
end

% Get new max voltage & check with tseries data vector that it's within limits.
tseries = getCurrentVoltage(handles);
if isnumeric( tseries.data )
   minv = min(tseries.data)*1.2;
   maxv = max(tseries.data)*1.2;
elseif iscell( tseries.data )
   minv = min( cellfun( @(d) min( d(:) ), tseries.data ) );
   maxv = max( cellfun( @(d) max( d(:) ), tseries.data ) );
end
str     = get(handles.voltage_max, 'String');
[ok, newv] = checkStringInput(str, 'float', minv, maxv);
if ~ok
   % users input is dodgy - gotta reverse engineer old text box max from
   % slider value & max display limit (zoom_max)
   % zoom_max   = text_max - (text_max - text_min)*proport;
   % zoom_max   = text_max - (text_max - text_min)*proport;
   percent = get(handles.voltage_slider, 'Value');
   proport = (1-percent)/2; % proportion max time changed by slider percentage
   zoom_max= handles.data.vlim(2);
   [~,text_min] = checkStringInput(get(handles.voltage_min, 'String'), 'float');
   text_max= (zoom_max + proport*text_min) / (1-proport);
   displayErrorMsg(sprintf('Voltage must be between %d & %d', minv, maxv));
   set(handles.voltage_max, 'String', sprintf('%.2f',text_max)); % reset to old val
   return
end
% Gotta get min time from other text box because need to reset displayed
% lims just in case new max is smaller than the last displayed min
% (zoom allows you to zoom in from text box mins/maxes)
[~,minv] = checkStringInput(get(handles.voltage_min, 'String'), 'float'); % we know str is valid
set(handles.voltage_max, 'String', sprintf('%.2f',newv)); % reset to new val
set(handles.voltage_slider, 'Value', 1);
handles.data.vlim(2) = newv;
handles.data.vlim(1) = minv;
handles = updateSETFigure(handles, tseries);
guidata(hObject, handles);
end

function voltage_min_Callback(hObject, eventdata, handles)
if ~haveUserData(handles)
   return;
end

% Get new max voltage & check with tseries data vector that it's within limits
tseries = getCurrentVoltage(handles);
minv    = min(tseries.data)*1.2;
maxv    = max(tseries.data)*1.2;
str     = get(handles.voltage_min, 'String');
[ok, newv] = checkStringInput(str, 'float', minv, maxv);
% if user's put in 0 assume they want the smallest voltage value, dt
if ~ok
   % users input is dodgy - gotta reverse engineer old text box min from
   % slider value & min display limit (zoom_min)
   % zoom_min = (text_max - text_min)*proport + text_min;
   % text_min = (zoom_min - text_max*proport) / (1-alpha)
   percent = get(handles.voltage_slider, 'Value');
   proport = (1-percent)/2; % proportion max voltage changed by slider percentage
   zoom_min= handles.data.vlim(1);
   [~,text_max] = checkStringInput(get(handles.voltage_max, 'String'), 'float');
   text_min= (zoom_min - proport*text_max) / (1-proport);
   displayErrorMsg(sprintf('Voltage must be between %d & %d', minv, maxv));
   set(handles.voltage_min, 'String', sprintf('%.2f',text_min)); % reset to old val
   return
end
% Gotta get max voltage from other text box because need to reset displayed
% voltage lims just in case new max is smaller than the last displayed min
% (zoom allows you to zoom in from text box mins/maxes)
[~,maxv] = checkStringInput(get(handles.voltage_max, 'String'), 'float'); % we know str is valid
set(handles.voltage_min, 'String', sprintf('%.2f',newv)); % reset to new val
set(handles.voltage_slider, 'Value', 1);
handles.data.vlim(1) = newv;
handles.data.vlim(2) = maxv;
handles = updateSETFigure(handles, tseries);
guidata(hObject, handles);
end

% --- Executes on button press in access_voltage.
function access_voltage_Callback(hObject, eventdata, handles)
[tseries, ~, ~, ts_name] = getCurrentVoltage(handles);
assignin('base', title2Str(ts_name,1,1,'_'), tseries);

end

% --- Executes on button press in new_figure.
function new_figure_Callback(hObject, eventdata, handles)
new_gui     = SpikeExtractionTool;
new_handles = guidata(gcf);
tseries     = getCurrentVoltage(handles);
handles.data.guihandles = [handles.data.guihandles; new_gui.guihandles]; % record new gui in our data struct
guidata(handles.figure1, handles);

%% Generate user data object to store
new_gui.num_tseries = 1;
new_gui.tseries     = {tseries};             % initialise cell array of time series
new_gui.curr_tseries= 1;                     % current time series
new_gui.tseries_str = {tseries.name};        % initialise names of time series
new_gui.used_names  = {{tseries.name}, 1};   % list of names already used by user
new_gui.last_tseries= 1;   % last time series - restore if current ts deleted
new_gui.tlim        = handles.data.tlim;
new_gui.vlim        = handles.data.vlim;
new_gui.last_tool   = 1;                     % record of user's last tool
new_gui.params      = getDefaultToolParams;  % default params for all tools & implementation methods
new_gui.last_dir    = handles.data.last_dir; % where they opened gui from
% copy of all SET gui handles - for each gui, make sure it's own handle is 1st
new_gui.guihandles  = [new_gui.guihandles; handles.data.guihandles(1:end-1)];
new_handles.data    = new_gui;               % record user data in handle
new_handles         = toggleSETGUIstate(new_handles, 'on'); % switch everything off until data's loaded

% udpate voltage timeseries drop down list
set(new_handles.curr_signal, 'String', new_handles.data.tseries_str);
set(new_handles.curr_signal, 'Value', 1);

% plot data that was put in new figure
curr_signal_Callback( new_handles.curr_signal, [], new_handles );
end

% --- Executes on slider movement.
function scroll_axes_Callback(hObject, eventdata, handles)
mouseWaitingFunction(handles.figure1, @scroll_axes, hObject, eventdata, handles);
end

function scroll_axes(hObject, eventdata, handles)
tseries = getCurrentVoltage(handles);
handles = updateSETFigure(handles, tseries);
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function scroll_axes_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% USER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    data.num_tseries = 0;
%    data.tseries     = [];   % cell array of tseries with name, dt, method used to generate etc
%    data.curr_tseries= [];   % current time series as an index into tseries cell array
%    data.tseries_str = [];   % initialise names of time series
%    data.used_names  = [];   % list of names already used by user
%    data.last_tseries= [];   % last time series - restore if current ts deleted
%    data.tlim        = [];   % record time limits chosen by user
%    data.vlim        = [];   % record time limits chosen by user
%    data.last_tool   = [];   % record of user's last tool
%    data.last_dir    = pwd;  % where they opened gui from
%    handles.data     = data; % record user data in handle

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

% Remove voltage from handles
function varargout = removeVoltage(handles, varargin)

if numel(varargin) == 0
   remove_tseries = handles.data.curr_tseries;
else
   remove_tseries = varargin{1}; % In case we are removing a voltage different to the currently active
end
   tseries_str    = handles.data.tseries_str{remove_tseries};
   last_tseries   = handles.data.last_tseries;
   num_tseries    = handles.data.num_tseries;

% if the last remaining timeseries is being deleted, clear gui
if num_tseries==1
   handles = toggleSETGUIstate(handles, 'off');
   guidata(handles.figure1, handles);
   if numel(varargin) == 1,varargout = {handles};end % When removing multiple voltages at once, we want to update the state of handles for the caller method
   return;
   
   % last time series may be the same as current time series (being removed)
   %  after deleting a timeseries or some such thing, so show something else
elseif handles.data.last_tseries==remove_tseries
   indices     = 1:num_tseries;
   indices(remove_tseries) = [];
   handles.data.last_tseries = indices(1);
end
% display last tseries before removing this one so we can check if data
% types are the same etc
set(handles.curr_signal, 'Value', handles.data.last_tseries);
% Update figure & handles structure
curr_signal_Callback(handles.curr_signal, [], handles);
handles = guidata(handles.figure1);
% curr_signal changes last & current timeseries indices, so change back
% (we wanted the function to display the previous timeseries, but the
% user didn't actually select it themselves so don't update current &
% last timeseries indices)
handles.data.last_tseries = last_tseries;
handles.data.curr_tseries = remove_tseries;

% now do the actual removing
handles.data.used_names = removeStringFromList(handles.data.used_names,  tseries_str);
handles.data.tseries(remove_tseries)     = [];
handles.data.tseries_str(remove_tseries) = [];
handles.data.num_tseries              = handles.data.num_tseries - 1;

if handles.data.num_tseries==0
   handles = toggleGUIstate(handles,'off');
   if numel(varargin) == 1,varargout = {handles};end % When removing multiple voltages at once, we want to update the state of handles for the caller method
   return
end
% update our indices for last tseries if necessary
if handles.data.last_tseries == remove_tseries
   handles.data.last_tseries = 1;
elseif handles.data.last_tseries > remove_tseries
   handles.data.last_tseries = handles.data.last_tseries - 1;
end
if handles.data.last_tseries > length(handles.data.tseries_str)
   str = sprintf('last timeseries index (%d) larger than number of strings (%d)',...
      handles.data.last_tseries, length(handles.data.tseries_str) );
   displayErrorMsg(str);
   handles.data.last_tseries = 1;
end
handles.data.curr_tseries = handles.data.last_tseries;
set(handles.curr_signal, 'String', handles.data.tseries_str);
set(handles.curr_signal, 'Value', handles.data.last_tseries);
handles.data.curr_tseries = handles.data.last_tseries;
guidata(handles.figure1, handles);

if numel(varargin) == 1,varargout = {handles};end % When removing multiple voltages at once, we want to update the state of handles for the caller method

end


% extract time series data, index number into tseries cell array, data type
% of time series, and the name of the time series
function [tseries, ts_num, type, ts_name] = getCurrentVoltage(handles)
% if try fails then we're doing if for viewVolume_newWindow - get stat
% from data struct instead
try
   ts_num  = get(handles.curr_signal,'Value');
   tseries = handles.data.tseries{ts_num};
   ts_name = handles.data.tseries_str{ts_num};
   type    = tseries.type;
catch
   ts_num  = handles.data.curr_tseries;
   tseries = handles.data.tseries{ts_num};
   ts_name = handles.data.tseries_str{ts_num};
   type    = tseries.type;
end
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

% function closeGUI(hObject, eventdata, handles) % is format below from old
% version of matlab??
function closeGUI(src, evnt)
diary off

handles  = guidata(src);
% if here we've got figure handles - check if there's user data, &
% quite immediately if there isn't
if isfield(handles, 'data')
   data = handles.data;
else
   return;
end
% don't get user confirmation if there's no user data!
if isfield(handles, 'data') && handles.data.num_tseries>0
   % get user confirmation that they want to close the gui
   response = userConfirmation('Do you want to close the GUI?',...
      'Close Request Function');
   switch response
      case 'Yes'
         % do nothing yet
      case 'No'
         diary on
         return;
   end
end

% if we've opened smr files we need to close the CED library
% unload the library when finished
try unloadlibrary ceds64int; catch, end


% NOTE: to get a list of all objects matching some criteria, but for
% guis, need to use findall() instead of findobj() if the handle
% visibility of the gui is set to 'callback'

% if we're here we're deleting the gui so need to remove it from the
% other gui's list of handles
try
   this_gui   = data.guihandles(1);
   other_guis = data.guihandles(2:end);
   for gi=1:length(other_guis)
      other = other_guis(gi);
      % make sure handle hasn't been deleted without us knowing somehow
      if ishandle(other)
         other_figures = guidata(other);
         % for some reason guidata can give either figure handle or user data
         if isfield(other_figures, 'data')
            other_data = other_figures.data;
         end
         other_handles = other_data.guihandles;
         if ~isempty(other_handles)
            deleted = ~ishandle(other_handles);
            other_handles(deleted) = [];
            tf = eq(this_gui, other_handles);
            other_handles(tf)  = [];
            other_data.guihandles = other_handles;
            guidata(other.figure1, other);
         end
      end
   end
catch ME
end

% no matter what happens with those goddam m.f. gui's, delete this one
delete(gcf);

end


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

% Mouse pointer back to normal.
set(handles.figure1, 'pointer', 'arrow')
end

%% Context menu options
% --------------------------------------------------------------------
function loadingWindow_Callback(hObject, eventdata, handles)
% hObject    handle to loadingWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if strcmp('on', hObject.Checked)
      hObject.Checked = 'off';
      handles.loadingWindowOn = false;
   else
      hObject.Checked = 'on';
      handles.loadingWindowOn = true;
   end
   % Update handles structure
   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function debugNone_Callback(hObject, eventdata, handles)
% hObject    handle to debugNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if strcmp('off', hObject.Checked)
      hObject.Checked = 'on';
      handles.debugFull.Checked = 'off';
      handles.debugFull.Enable = true;
      handles.debugSemi.Checked = 'off';
      handles.debugSemi.Enable = true;
      handles.debugOption = 'none';
      hObject.Enable = false;
   end
   % Update handles structure
   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function debugSemi_Callback(hObject, eventdata, handles)
% hObject    handle to debugSemi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if strcmp('off', hObject.Checked)
      hObject.Checked = 'on';
      handles.debugFull.Checked = 'off';
      handles.debugFull.Enable = true;
      handles.debugNone.Checked = 'off';
      handles.debugNone.Enable = true;
      handles.debugOption = 'semi';
      hObject.Enable = false;
   end
   % Update handles structure
   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function debugFull_Callback(hObject, eventdata, handles)
% hObject    handle to debugFull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if strcmp('off', hObject.Checked)
      hObject.Checked = 'on';
      handles.debugSemi.Checked = 'off';
      handles.debugSemi.Enable = true;
      handles.debugNone.Checked = 'off';
      handles.debugNone.Enable = true;
      handles.debugOption = 'full';
      hObject.Enable = false;
   end   
   % Update handles structure
   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function rescaleAtTheEnd_Callback(hObject, eventdata, handles)
% hObject    handle to rescaleAtTheEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if strcmp('off', hObject.Checked)
      hObject.Checked = 'on';
      handles.rescalePeaks.Checked = 'off';
      handles.rescalePeaks.Enable = true;
      handles.rescaleJumpAhead.Checked = 'off';
      handles.rescaleJumpAhead.Enable = true;
      handles.rescaleOption = 'at_end';
      hObject.Enable = false;
   end   
   % Update handles structure
   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function rescalePeaks_Callback(hObject, eventdata, handles)
% hObject    handle to rescalePeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if strcmp('off', hObject.Checked)
      hObject.Checked = 'on';
      handles.rescaleAtTheEnd.Checked = 'off';
      handles.rescaleAtTheEnd.Enable = true;
      handles.rescaleJumpAhead.Checked = 'off';
      handles.rescaleJumpAhead.Enable = true;
      handles.rescaleOption = 'peaks';
      hObject.Enable = false;
   end   
   % Update handles structure
   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function rescaleJumpAhead_Callback(hObject, eventdata, handles)
% hObject    handle to rescaleJumpAhead (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if strcmp('off', hObject.Checked)
      hObject.Checked = 'on';
      handles.rescaleAtTheEnd.Checked = 'off';
      handles.rescaleAtTheEnd.Enable = true;
      handles.rescalePeaks.Checked = 'off';
      handles.rescalePeaks.Enable = true;
      handles.rescaleOption = 'JA';
      hObject.Enable = false;
   end   
   % Update handles structure
   guidata(hObject, handles);
end


% --------------------------------------------------------------------
function clearDifferentMenu_Callback(hObject, eventdata, handles)
% hObject    handle to clearDifferentMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if handles.data.num_tseries == 1
      displayErrorMsg('There is only one voltage loaded.');
      return
   else
      % Find position of current mouse
      get(handles.figure1,'Parent');
      screenHandle = get(handles.figure1,'Parent');
      pos = screenHandle.PointerLocation;
      mainData = guidata(gcbf);
      clearVoltageWindow = figure('Name','Clear','MenuBar','none',...
         'ToolBar','none','Position',[pos(1),pos(2)-266,212,266],...
         'WindowStyle','modal','Resize','off','DockControls','off');
      vtc = uicontrol('parent',clearVoltageWindow,'Style','listbox',...
         'String',mainData.data.tseries_str,'InnerPosition',[14,56,187,155],...
         'Enable','on','Max',mainData.data.num_tseries,'Tag','vtc');
      uicontrol('parent',clearVoltageWindow,'Style','text',...
         'String', 'Select the voltage(s) you want to remove:',...
         'InnerPosition',[14,211,187,50],'FontSize',11);
      uicontrol('parent',clearVoltageWindow,'Style','pushbutton',...
         'String','Clear','Position',[130,15,69,30],'Callback',...
         @(hObject,eventdata)SpikeExtractionTool('clearVoltages',hObject,eventdata,mainData,vtc));
   end
end

function clearVoltages(hObject, eventdata, mainData, vtc)
   % Set the mouse pointer to waiting to know the function is running.
   set(mainData.figure1, 'pointer', 'watch')
   drawnow;
   try
      vtc = vtc.Value;
      delete(hObject.Parent);
      for i = numel(vtc):-1:1
         mainData = removeVoltage(mainData, vtc(i));
      end
   catch E
      % Mouse pointer back to normal.
      set(mainData.figure1, 'pointer', 'arrow')
      runtimeErrorHandler(E);
   end
   % Mouse pointer back to normal.
   set(mainData.figure1, 'pointer', 'arrow')
end


% --------------------------------------------------------------------
function exitMenuBar_Callback(hObject, eventdata, handles)
% hObject    handle to exitMenuBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
closeGUI(hObject,eventdata);
end


% --------------------------------------------------------------------
function plotNewFigMenuBar_Callback(hObject, eventdata, handles)
% hObject    handle to plotNewFigMenuBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_figure_Callback(hObject,eventdata,handles);
end


% --------------------------------------------------------------------
function accessVoltageMenuBar_Callback(hObject, eventdata, handles)
% hObject    handle to accessVoltageMenuBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
access_voltage_Callback(hObject, eventdata, handles);
end


% --------------------------------------------------------------------
function saveVoltageMenuBar_Callback(hObject, eventdata, handles)
% hObject    handle to saveVoltageMenuBar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_voltage_Callback(hObject, eventdata, handles);
end

function removeToolBarButtons()
removeItems = ([{'Save Figure'}, {'New Figure'}, {'Open File'}, {'Print Figure'}, {'Link Plot'}, {'Open Property Inspector'}, {'Insert Colorbar'}]);
for i = 1:size(removeItems,2)
   listOfElements = findall(gcf);
   element = findall(listOfElements,'ToolTipString',string(removeItems(i)));
   set(element,'Visible','off');
end
end


% --- Executes on button press in toggleZoomButton.
function toggleZoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to toggleZoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Toggles between arrows (displacement) and magnifier (zoom) icons
   warning('off','MATLAB:imagesci:png:libraryWarning'); % Ignore PNG associated warning

   if strcmp('zoom',hObject.UserData)
      % Displacement function has been selected
      [x,~]=imread('fig/arrowsIcon.png');% Load the displacement icon
      I2=imresize(x, [22 22]); % Resize icon
      hObject.CData = I2; % Assign icon to the button
      hObject.UserData = 'disp'; % Change state to displacement
      handles.time_slider.Value = handles.data.displacementPercentage(1); % Update the position of the slider to represent displacement.
      handles.voltage_slider.Value = handles.data.displacementPercentage(2); % Update the position of the slider to represent displacement.
      handles.time_slider.SliderStep(2) = max([handles.time_slider.SliderStep(1) 1e-1/handles.data.zoomPercentage(1)]); % Change the size of the horizontal slider indicator to match the value zoomed in.
      handles.time_slider.SliderStep(2) = min(1, handles.time_slider.SliderStep(2)); % Make sure it is within 0 and 1.
      handles.time_slider.SliderStep(1) = (0.1 * handles.time_slider.SliderStep(2));
      handles.voltage_slider.SliderStep(2) = max(handles.voltage_slider.SliderStep(1) , handles.data.zoomPercentage(2));% Change the size of the vertical slider indicator to match the value zoomed in.
      handles.voltage_slider.SliderStep(1) = 0.1 * handles.voltage_slider.SliderStep(2);
   else
      % Zoom function has been selected
      [x,~]=imread('fig/magnifierIcon.png');% Load the zoom icon
      I2=imresize(x, [22 22]); % Resize icon
      hObject.CData = I2; % Assign icon to the button
      hObject.UserData = 'zoom'; % Change state to zoom
      handles.time_slider.Value = handles.data.zoomPercentage(1); % Update the position of the slider to represent zoom.
      handles.voltage_slider.Value = handles.data.zoomPercentage(2); % Update the position of the slider to represent zoom.
      handles.time_slider.SliderStep =  [0.001 0.1];%max(handles.time_slider.SliderStep(1) , handles.data.displacementPercentage(1)); % Change the size of the horizontal slider indicator to match the value displaced in.
      handles.voltage_slider.SliderStep =  [0.01 0.1];%max(handles.voltage_slider.SliderStep(1) , handles.data.displacementPercentage(2)); %  Change the size of the vertical slider indicator to match the value displaced in.
   end
end


% --- Executes during object creation, after setting all properties.
function toggleZoomButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to toggleZoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
   warning('off','MATLAB:imagesci:png:libraryWarning'); % Ignore PNG associated warning
   [x,map]=imread('fig/magnifierIcon.png'); % Load the zoom icon
   I2=imresize(x, [22 22]); % Resize icon
   hObject.CData = I2; % Assign icon to the button
   hObject.BackgroundColor = [1 0.6 0.6]; % Change color to match other buttons
   hObject.String = ''; % Remove string
   hObject.UserData = 'zoom'; % Change state to zoom     
end
