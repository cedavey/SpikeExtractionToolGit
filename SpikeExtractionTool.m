%% Copyright 2019 Artemio Soto-Breceda
% 
% Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3.0 ("the License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     https://www.gnu.org/licenses/gpl-3.0-standalone.html
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

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
% - also when extracting spikes a one-off large spike spawns a new family,
%   which it perhaps shouldn't? Perhaps make allowed size diff. dependent
%   on noise?
% - button to reset axes (& perhaps to refresh?)
% - implement mean-shift clustering & PCA to separate neurons
% - if I find crazy voltage values, should I set them to 0?
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
%    load_voltage_button    - press button to load new data file & erase current data
%    add_voltage_button     - button to add a new data file but retain all current data
%  Display voltage
%    curr_signal     - pop up list to choose current time series displayed
%    clear_voltage_button   - button to delete current voltage time series, but keep the rest
%  Generation parameters
%    param1(2,3...)  - text display of name of parameter used to generate display tseries
%    val1(2,3,...)   - text display of param value used to generate display tseries
%  Available tools
%    select_tool     - text display of select tool instruction
%    tool_list       - popup menu showing tools to apply to displayed tseries
%    select_method   - text display of select method instruction
%    method_list     - popup menu showing methods to apply current tool to displayed tseries
%    run_tool_button        - initiate implementation of tool using select method
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
%    scroll_axes_slider     - if have more axes in figure than can be viewed allow scrolling
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

% Last Modified by GUIDE v2.5 31-Jan-2022 12:09:02

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

% Create an instance of frontEndFunctions (shared functions between GUIDE
% and App) to access the functions.
handles.f = frontEndFunctions('gui');

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

handles          = toggleSETGUIstate(handles,'off'); % switch everything off until data's loaded

% assign gui a name
set(hObject, 'Name', 'Spike Extraction Tool');
% set(gcf, 'menubar', 'figure') % turn the figure menu on (so can save etc) % I don't think it is a good idea to have this on, better to make a custom menu bar.


% prepare for user confirmation when GUI is closed
set(handles.figure1, 'CloseRequestFcn', @closeGUI);
%     set(gcf, 'WindowScrollWheelFcn', {@figure1_figScroll, handles});
removeToolBarButtons(handles);

% from property inspector in guide, under WindowScrollWheelFcn
% @(hObject,eventdata)SpikeExtractionTool('figure1_WindowScrollWheelFcn',hObject,eventdata,guidata(hObject))

% Initialize options for debug and loading window
handles.options.loadingWindowOn = true;
handles.options.mouseWaitingOn = true;
handles.options.debugOption = 'none';
handles.options.rescaleOption = 'at_end'; % Within one jump ahead, what
	% segment of the timeseries is rescaled, 'at_end' (default) = rescles
    % after processing whole time series, and the rescaling periods are 
    % the "jump ahead" segments, 'peaks' = rescales in real time, from one
    % spike peak to the next one,'JA' = rescales in real time within each 
    % "jump ahead" segment.
handles.options.auto_params = 'false';

% Initialize tooltips
handles = setTooltips(handles);

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

% Get location of log files
path = getFilePath('log');
fid  = fopen([path 'log_all.log'], 'a'); % Opens log file to append this session's string
fprintf(fid, '\n\n-------------- %s @ %s | %s ---------------\n', getenv('Username'),getenv('UserDomain'),datestr(now, 0));
fclose(fid); % Close log file
diary([path 'log_all.log']); % Activates the diary function, i.e. save all the activity into a file.
handles.options.logPath = path;

% Update handles structure
guidata(hObject, handles);
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
if ~handles.f.haveUserData(handles), return; end

%     % see if pointer's in one of the axes
%     in_cp = isPointerInObject(handles, 'axes_coronal', 'crosshair');
%     if in_cp, return; end
%
%     % we're nowhere special so display pointer as an arrow
%     set(hObject,'pointer','arrow');
end

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
if ~handles.f.haveUserData(handles), return; end

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
tseries = handles.f.getCurrentVoltage(handles);
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

% --- Executes on button press in load_voltage_button.
function load_voltage_button_Callback(hObject, eventdata, handles)
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


if handles.options.mouseWaitingOn
    mouseWaitingFunction(handles.figure1,@load_voltage,hObject,eventdata,handles);
else
    load_voltage(hObject,eventdata,handles);
end
end

function load_voltage(hObject, eventdata, handles)
   handles.f.load_voltage(hObject, handles);
end

% --- Executes on selection change in curr_signal.
function curr_signal_Callback(hObject, eventdata, handles)
    if handles.options.mouseWaitingOn
        mouseWaitingFunction(handles.figure1,@different_signal_selected,hObject,eventdata,handles);
    else
        different_signal_selected(hObject,eventdata,handles);
    end
end

function different_signal_selected(hObject, eventdata, handles)
   handles.f.curr_signal(hObject,eventdata,handles);
end

% --- Executes during object creation, after setting all properties.
function curr_signal_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', '');
set(hObject, 'Value', 1);
end

% --- Executes on button press in add_voltage_button.
function add_voltage_button_Callback(hObject, eventdata, handles)
    if handles.options.mouseWaitingOn
        mouseWaitingFunction(handles.figure1,@add_voltage,hObject,eventdata,handles);
    else
        add_voltage(hObject,eventdata,handles);
    end
end

function add_voltage(hObject,eventdata,handles)
   handles.f.add_voltage(handles);
end

% --- Executes on button press in clear_voltage_button.
function clear_voltage_button_Callback(hObject, eventdata, handles)
   [tseries, ts_num, data_type, ts_name] = handles.f.getCurrentVoltage(handles);
   response = userConfirmation(['Delete ' ts_name '?'],...
      'Clear current time series?');
   if strcmp(response,'No')
      return
   end
   
   if handles.options.mouseWaitingOn
        mouseWaitingFunction(handles.figure1,@removeVoltage,handles); % Instead of removeVoltage(handles);
    else
        removeVoltage(handles);
    end
end

% --- Executes on selection change in tool_list.
function tool_list_Callback(hObject, eventdata, handles)
   handles.f.tool_list(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function tool_list_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in method_list.
function method_list_Callback(hObject, eventdata, handles)
   % update the tooltips
   tooltip_name = {['method_list_' lower(handles.method_list.String{handles.method_list.Value})]};
   handles = setTooltips(handles, {'method_list'}, getTooltips(tooltip_name));
end

% --- Executes during object creation, after setting all properties.
function method_list_CreateFcn(hObject, eventdata, handles)
   if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
      set(hObject,'BackgroundColor','white');
   end
end

% --- Executes on button press in set_tool_params.
function set_tool_params_Callback(hObject, eventdata, handles)
   handles.f.set_tool_params(handles, hObject);
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

% if removing templates. User has chosen how many to remove.
function userids = getUserTemplateDeleteIDs(handles, tseries, method_params )
   userids = handles.f.getUserTemplateDeleteIDs(tseries, method_params);
end

% --- Executes on button press in run_tool_button.
function run_tool_button_Callback(hObject, eventdata, handles)
    if handles.options.mouseWaitingOn
        mouseWaitingFunction(handles.figure1,@run_tool,hObject,eventdata,handles);
    else
        run_tool(hObject,eventdata,handles);
    end
end

function run_tool(hObject, eventdata, handles)
   handles.f.run_tool(handles);
end

% --- Executes on button press in save_voltage.
function save_voltage_Callback(hObject, eventdata, handles)
   handles.f.save_voltage(handles, hObject);
end

% --- Executes on slider movement.
function time_slider_Callback(hObject, eventdata, handles)
    if ~handles.f.haveUserData(handles)
    	return;
    end
   
    if handles.options.mouseWaitingOn
        mouseWaitingFunction(handles.figure1,@time_slider_updated,hObject,eventdata,handles);
    else
        time_slider_updated(hObject,eventdata,handles);
    end
end

function time_slider_updated(hObject,eventdata,handles)
   handles.f.time_slider_updated(handles, hObject);
end


% --- Executes during object creation, after setting all properties.
function time_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function time_max_Callback(hObject, eventdata, handles)
   handles.f.time_max(handles,hObject);
end

function time_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
end

function time_min_Callback(hObject, eventdata, handles)
   handles.f.time_min(handles, hObject);
end


% --- Executes during object creation, after setting all properties.
function time_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor','white');
end
end

% --- Executes on slider movement.
function voltage_slider_Callback(hObject, eventdata, handles)
   if ~handles.f.haveUserData(handles)
      return;
   end

   if handles.options.mouseWaitingOn
       mouseWaitingFunction(handles.figure1, @voltage_slider_updated, handles, hObject);
   else
       voltage_slider_updated(handles, hObject);
   end
end

function voltage_slider_updated(handles, hObject)
   handles.f.voltage_slider_updated(handles, hObject);
end

% --- Executes during object creation, after setting all properties.
function voltage_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
   set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function voltage_max_Callback(hObject, eventdata, handles)

   handles.f.voltage_max(handles, hObject);

% if ~handles.f.haveUserData(handles)
%    return;
% end
% 
% % Get new max voltage & check with tseries data vector that it's within limits.
% tseries = handles.f.getCurrentVoltage(handles);
% if isnumeric( tseries.data )
%    minv = min(tseries.data)*1.2;
%    maxv = max(tseries.data)*1.2;
% elseif iscell( tseries.data )
%    try
%       minv = min( cellfun( @(d) min( d(:) ), tseries.data ) );
%       maxv = max( cellfun( @(d) max( d(:) ), tseries.data ) );
%    catch
%       % If there is an error, it means the cell structure is a bit more
%       % complicated. Set the min and max based on the first template
%       minv = min(min(tseries.data{1}{1}.spikes)) * 1.2;
%       maxv = max(max(tseries.data{1}{1}.spikes)) * 1.2;
%    end
% end
% str     = get(handles.voltage_max, 'String');
% [ok, newv] = checkStringInput(str, 'float', minv, maxv);
% if ~ok
%    % users input is dodgy - gotta reverse engineer old text box max from
%    % slider value & max display limit (zoom_max)
%    % zoom_max   = text_max - (text_max - text_min)*proport;
%    % zoom_max   = text_max - (text_max - text_min)*proport;
%    percent = get(handles.voltage_slider, 'Value');
%    proport = (1-percent)/2; % proportion max time changed by slider percentage
%    zoom_max= handles.data.vlim(2);
%    [~,text_min] = checkStringInput(get(handles.voltage_min, 'String'), 'float');
%    text_max= (zoom_max + proport*text_min) / (1-proport);
%    displayErrorMsg(sprintf('Voltage must be between %d & %d', minv, maxv));
%    set(handles.voltage_max, 'String', sprintf('%.2f',text_max)); % reset to old val
%    return
% end
% % Gotta get min time from other text box because need to reset displayed
% % lims just in case new max is smaller than the last displayed min
% % (zoom allows you to zoom in from text box mins/maxes)
% [~,minv] = checkStringInput(get(handles.voltage_min, 'String'), 'float'); % we know str is valid
% set(handles.voltage_max, 'String', sprintf('%.2f',newv)); % reset to new val
% 
% handles.data.vlim(2) = newv;
% if ~(minv < newv)
%    currtseries = handles.data.curr_tseries;
%    minv = min(handles.data.tseries{currtseries}.data);
%    handles.voltage_min.String = num2str(minv);
% end
% handles.data.vlim(1) = minv;
% 
% % Update voltage sliders
% handles.f.updateVoltageSlider(handles);
% 
% handles = updateSETFigure(handles, tseries);
% guidata(hObject, handles);
end

function voltage_min_Callback(hObject, eventdata, handles)
   handles.f.voltage_min(handles, hObject);
end

% --- Executes on button press in access_voltage.
function access_voltage_Callback(hObject, eventdata, handles)
[tseries, ~, ~, ts_name] = handles.f.getCurrentVoltage(handles);
assignin('base', title2Str(ts_name,1,1,'_'), tseries);

end

% --- Executes on button press in new_figure.
function new_figure_Callback(hObject, eventdata, handles)
new_gui     = SpikeExtractionTool;
new_handles = guidata(gcf);
tseries     = handles.f.getCurrentVoltage(handles);
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
% curr_signal_changed( new_handles.curr_signal, [], new_handles );
new_handles.f.curr_signal(new_handles.curr_signal, '', handles);
new_handles = updateSETFigure(new_handles, tseries);
guidata(new_handles.figure1, new_handles);
end

% --- Executes on slider movement.
function scroll_axes_slider_Callback(hObject, eventdata, handles)
if handles.options.mouseWaitingOn
    mouseWaitingFunction(handles.figure1, @scroll_axes, hObject, eventdata, handles);
else
    scroll_axes(hObject,eventdata,handles);
end
end

function scroll_axes(hObject, eventdata, handles)
tseries = handles.f.getCurrentVoltage(handles);
handles = updateSETFigure(handles, tseries);
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function scroll_axes_slider_CreateFcn(hObject, eventdata, handles)
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

function name = getFileName(handles,instruct, defVal, maxLength)
   name = handles.f.getFileName(instruct, defVal, maxLength);
end

% Remove voltage from handles

function varargout = removeVoltage(handles, varargin)
   try
      varargout = handles.f.removeVoltage(handles, varargin);
   catch E
      if strcmp('MATLAB:unassignedOutputs', E.identifier)
         varargout = {};
      else
         rethrow(E);
      end
   end
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
      handles.options.loadingWindowOn = false;
   else
      hObject.Checked = 'on';
      handles.options.loadingWindowOn = true;
   end
   % Update handles structure
   guidata(hObject, handles);
end

% --------------------------------------------------------------------
function mouseWaiting_Callback(hObject, eventdata, handles)
   if strcmp('on', hObject.Checked)
      hObject.Checked = 'off';
      handles.options.mouseWaitingOn = false;
   else
      hObject.Checked = 'on';
      handles.options.mouseWaitingOn = true;
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
      handles.options.debugOption = 'none';
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
      handles.options.debugOption = 'semi';
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
      handles.options.debugOption = 'full';
      hObject.Enable = false;
   end
   % Update handles structure
   guidata(hObject, handles);
end
% 
% % --------------------------------------------------------------------
% function rescaleAtTheEnd_Callback(hObject, eventdata, handles)
% % hObject    handle to rescaleAtTheEnd (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%    if strcmp('off', hObject.Checked)
%       hObject.Checked = 'on';
%       handles.rescalePeaks.Checked = 'off';
%       handles.rescalePeaks.Enable = true;
%       handles.rescaleJumpAhead.Checked = 'off';
%       handles.rescaleJumpAhead.Enable = true;
%       handles.options.rescaleOption = 'at_end';
%       hObject.Enable = false;
%    end
%    % Update handles structure
%    guidata(hObject, handles);
% end
% 
% % --------------------------------------------------------------------
% function rescalePeaks_Callback(hObject, eventdata, handles)
% % hObject    handle to rescalePeaks (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%    if strcmp('off', hObject.Checked)
%       hObject.Checked = 'on';
%       handles.rescaleAtTheEnd.Checked = 'off';
%       handles.rescaleAtTheEnd.Enable = true;
%       handles.rescaleJumpAhead.Checked = 'off';
%       handles.rescaleJumpAhead.Enable = true;
%       handles.options.rescaleOption = 'peaks';
%       hObject.Enable = false;
%    end
%    % Update handles structure
%    guidata(hObject, handles);
% end
% 
% % --------------------------------------------------------------------
% function rescaleJumpAhead_Callback(hObject, eventdata, handles)
% % hObject    handle to rescaleJumpAhead (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%    if strcmp('off', hObject.Checked)
%       hObject.Checked = 'on';
%       handles.rescaleAtTheEnd.Checked = 'off';
%       handles.rescaleAtTheEnd.Enable = true;
%       handles.rescalePeaks.Checked = 'off';
%       handles.rescalePeaks.Enable = true;
%       handles.options.rescaleOption = 'JA';
%       hObject.Enable = false;
%    end
%    % Update handles structure
%    guidata(hObject, handles);
% end


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
   mainData.f.clear_voltages(hObject, mainData, vtc);
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

function removeToolBarButtons(handles)

   v = ver('matlab');
   % Check if version is older than 2018b
   if str2double(v.Version) >= 9.5
     removeToolbarExplorationButtons(handles.figure1.Children(12));
   end
   removeItems = ([{'Save Figure'}, {'Insert Legend'}, {'Edit Plot'}, {'New Figure'}, {'Open File'}, {'Print Figure'}, {'Link Plot'}, {'Open Property Inspector'}, {'Insert Colorbar'}]);
   for i = 1:size(removeItems,2)
      listOfElements = findall(gcf);
      element = findall(listOfElements,'ToolTipString',string(removeItems(i)));
      set(element,'Visible','off');
   end
end

% --- Executes on button press in toggleZoomButton.
function toggleZoomButton_Callback(hObject, eventdata, handles)
   handles.f.toggleZoomButton(handles, hObject);
end

% --- Executes during object creation, after setting all properties.
function toggleZoomButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to toggleZoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
   warning('off','MATLAB:imagesci:png:libraryWarning'); % Ignore PNG associated warning

   % Get location of log files
   path    = getFilePath();
   [x,map] = imread([path 'fig' filesep 'magnifierIcon.png']); % Load the zoom icon
   I2      = imresize(x, [22 22]); % Resize icon
   
   hObject.CData = I2; % Assign icon to the button
   hObject.BackgroundColor = [1 0.6 0.6]; % Change color to match other buttons
   hObject.String = ''; % Remove string
   hObject.UserData = 'zoom'; % Change state to zoom
end


% --------------------------------------------------------------------
function helpMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to helpMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   % Put the file on the web so don't have to recompile whole gui to change help file
   %path = getFilePath();
   %open_pdf([path 'resources' filesep 'SEThelp.pdf']);
   web('https://github.com/srarty/SpikeExtractionToolGit/files/3823676/SEThelp.pdf','-browser');
end

% --------------------------------------------------------------------
function aboutMenuItem_Callback(hObject, eventdata, handles)
   handles.f.aboutMenuItem(handles);
end


% --- Executes on button press in automatic_params.
function automatic_params_Callback(hObject, eventdata, handles)
   handles.f.automatic_params(handles, hObject);
end


% --------------------------------------------------------------------
function batchProcessingMenu_Callback(hObject, eventdata, handles)
   try
      if handles.options.mouseWaitingOn
        mouseWaitingFunction(handles.figure1,@runBatchProcessing,hObject,eventdata,handles);
      else
        runBatchProcessing(hObject,eventdata,handles);
      end
   catch E
      handles.options.isBatch = false;
      guidata(hObject,handles);
      str = sprintf('\tUnexpected error during batch processing\n');
      caught_error = runtimeErrorHandler(E, 'message', str);
      if ~isempty(caught_error), rethrow(caught_error); end
   end
end

function runBatchProcessing(hObject, eventdata, handles)
% hObject    handle to batchProcessingMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   handles.f.runBatchProcessing(handles, hObject, eventdata);
end

function figure1_ButtonDownFcn(hObject, eventdata, handles)

end


% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles)
    if handles.options.mouseWaitingOn
        mouseWaitingFunction(handles.figure1, @reset_button_clicked, handles, hObject);
    else
        reset_button_clicked(handles, hObject);
    end
end

function reset_button_clicked(handles, hObject)
   handles.f.resetButtonClicked(handles, hObject);
end

% --- Executes during object creation, after setting all properties.
function reset_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
   warning('off','MATLAB:imagesci:png:libraryWarning'); % Ignore PNG associated warning

   % Get location of files
   % Get location of log files
   path = getFilePath();

   [x,map]=imread([path 'fig' filesep 'resetIcon.png']); % Load the zoom icon
   I2=imresize(x, [22 22]); % Resize icon
   hObject.CData = I2; % Assign icon to the button
   hObject.BackgroundColor = [1 0.6 0.6]; % Change color to match other buttons
   hObject.String = ''; % Remove string
end


% --- Executes on button press in new_figure_workspace.
% Loads the voltage to a workspace variable and plots it, then clears
% the variable
function new_figure_workspace_Callback(hObject, eventdata, handles)
    % Load the voltage to the workspace
    try
        [tseries, ~, ~, ts_name] = handles.f.getCurrentVoltage(handles);
        assignin('base', title2Str(ts_name,1,1,'_'), tseries);
        % Plot the figure:
        switch tseries.type
            case 'voltage'
                h = figure('Name', ts_name);
                plot(tseries.time, tseries.data);
                xlabel('Time (s)');
                ylabel('Vm (mV)');
                box off;
            case 'spike'
                h = figure('Name', ts_name);
                plot(tseries.time, tseries.data{1});
                if numel(tseries.data) > 1
                    hold on;
                    for i = 2:numel(tseries.data)
                        plot(tseries.time, tseries.data{i});
                    end
                end                        
                xlabel('Time (s)');
                ylabel('Vm (mV)');
                box off;
            otherwise
                error('Not implemented');
        end
    catch ME
        if exist('h', 'var')
            close(h);
            clear h
        end
       str = getCatchMEstring( ME, 'Error trying to plot the voltage.', false );
       displayErrorMsg( 'Sorry! I don''t know how to plot that kind of data yet. :(');
    end
end




% --------------------------------------------------------------------
function simulator_button_Callback(hObject, eventdata, handles)
% hObject    handle to simulator_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    try
        str = 'Looking for SpikeSimulationTool in the current parent folder.\n';
        printMessage('off','Text',str);
        d = dir(['..' filesep]);
        for i = 1:size(d)
            if strcmp('SpikeSimulationTool1',d(i).name)
                str = 'Attempting to run SST.m\n';
                printMessage('off','Text',str);
                addpath(['..' filesep 'SpikeSimulationTool']);
                run('SST');
                str = 'Done.\n';
                printMessage('off','Keyword',str);
                return
            end
        end
        % If reaches this point, it didn't find the Simulator folder with
        % that name. Allow user to find it manually.
        str = 'Folder ''SpikeSimulationTool'' not found. Manual entry to find the path.\n';
        printMessage('off','Text',str);
        d = uigetdir(['..' filesep]);
        addpath(d);
        run('SST');
    catch E
        error(E);
    end
end
