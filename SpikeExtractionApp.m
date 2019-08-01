classdef SpikeExtractionApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        f                      frontEndFunctions
        figure1                matlab.ui.Figure
        fileMenu               matlab.ui.container.Menu
        loadMenuBar            matlab.ui.container.Menu
        addVoltageMenuBar      matlab.ui.container.Menu
        clearMenuBar           matlab.ui.container.Menu
        clearSelectedMenuBar   matlab.ui.container.Menu
        clearDifferentMenuBar  matlab.ui.container.Menu
        plotNewFigMenuBar      matlab.ui.container.Menu
        accessVoltageMenuBar   matlab.ui.container.Menu
        saveVoltageMenuBar     matlab.ui.container.Menu
        exitMenuBar            matlab.ui.container.Menu
        helpMenu               matlab.ui.container.Menu
        helpMenuItem           matlab.ui.container.Menu
        aboutMenuItem          matlab.ui.container.Menu
        load_panel             matlab.ui.container.Panel
        load_voltage_button    matlab.ui.control.Button
        select_dataset         matlab.ui.control.Label
        add_voltage_button     matlab.ui.control.Button
        voltage_panel          matlab.ui.container.Panel
        display_voltage        matlab.ui.control.Label
        curr_signal            matlab.ui.control.DropDown
        clear_voltage_button   matlab.ui.control.Button
        tool_panel             matlab.ui.container.Panel
        tools_title            matlab.ui.control.Label
        tool_list              matlab.ui.control.DropDown
        select_tool            matlab.ui.control.Label
        method_list            matlab.ui.control.DropDown
        select_method          matlab.ui.control.Label
        set_tool_params        matlab.ui.control.Button
        run_tool_button        matlab.ui.control.Button
        automatic_params       matlab.ui.control.CheckBox
        param_panel            matlab.ui.container.Panel
        gen_params             matlab.ui.control.Label
        param1                 matlab.ui.control.Label
        val1                   matlab.ui.control.Label
        param2                 matlab.ui.control.Label
        val2                   matlab.ui.control.Label
        param3                 matlab.ui.control.Label
        val3                   matlab.ui.control.Label
        param4                 matlab.ui.control.Label
        val4                   matlab.ui.control.Label
        param5                 matlab.ui.control.Label
        val5                   matlab.ui.control.Label
        param6                 matlab.ui.control.Label
        val6                   matlab.ui.control.Label
        param7                 matlab.ui.control.Label
        val7                   matlab.ui.control.Label
        param8                 matlab.ui.control.Label
        val8                   matlab.ui.control.Label
        save_voltage           matlab.ui.control.Button
        plot_panel             matlab.ui.container.Panel
        axes1                  matlab.ui.control.UIAxes
        voltage_slider         matlab.ui.control.Slider
        voltage_max            matlab.ui.control.EditField
        voltage_min            matlab.ui.control.EditField
        time_slider            matlab.ui.control.Slider
        time_max               matlab.ui.control.EditField
        time_min               matlab.ui.control.EditField
        access_voltage         matlab.ui.control.Button
        new_figure             matlab.ui.control.Button
        scroll_axes_slider     matlab.ui.control.Slider
        toggleZoomButton       matlab.ui.control.Button
        data
        options
    end
    
    methods (Access = private)
        function add_voltage(app)
            app.f.add_voltage(app);
        end
        
        function clearVoltages(app, hObject, event, mainData, vtc)
           app.f.clear_voltages(event.Source, mainData, vtc);
        end
        
        function closeGUI(src, eventdata)
            diary off
            
            handles  = guidata(src);
            % if here we've got figure handles - check if there's user data, &
            % quite immediately if there isn't
            if isfield(handles, 'data')
               data = app.data;
            else
               return;
            end
            % don't get user confirmation if there's no user data!
            if isfield(handles, 'data') && app.data.num_tseries>0
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
        
        function figure1_CloseRequestFcn(app)
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
            if app.data.num_tseries == 0
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
        
        function figure1_WindowButtonDownFcn(app, hObject, eventdata, handles)
            if ~haveUserData(handles), return; end
        end
        
        function figure1_WindowButtonMotionFcn(app, hObject, eventdata, handles)
            if ~haveUserData(handles), return; end
        end
                
        function name = getFileName(app, instruct, defVal, maxLength)
            name = app.f.getFileName(instruct, defVal, maxLength);
        end
        
        function userids = getUserTemplateDeleteIDs(app, tseries, method_params)
           % TODO - validate method works
            userids = app.f.getUserTemplateDeleteIDs(tseries, method_params);
        end
        
        function userids = getUserTemplateMergeIDs(app, tseries, method_params)
           % TODO - validate method works
            userids = app.f.getUserTemplateMergeIDs(tseries, method_params);
        end
       
        function inobject = isPointerInObject(app, handles, object, pointer)
            % E.g. position = get(app.voxelValue_top, 'Position');
            estr = ['pos = get(app.' object ', ''Position'');'];
            eval(estr);
            pt   = get(app.figure1, 'CurrentPoint'); % current pointer position
            
            inobject = (pt(1,1)>=pos(1) && pt(1,1)<=(pos(1)+pos(3))) && ...
               (pt(1,2)>=pos(2) && pt(1,2)<=(pos(2)+pos(4)));
            if inobject && exist('pointer', 'var')
               set(app.figure1,'pointer',pointer);
            end
        end
        
        function load_voltage(app)
            app.f.load_voltage(app);
        end
        
        function handles = move_crosshairs(app, handles, axis)
            % %     data_struct = app.data_struct; % mem probs - avoid copying
            tseries = app.f.getCurrentVoltage(app);
            maxt    = get(app.time_slider,   'Max');
            maxv    = get(app.voltage_slider,'Max');
            
            switch axis
               case 'lateral'
                  %             % get new crosshair position values
                  %             buttonpnt = get(app.axes_lateral,'CurrentPoint');
                  %             slice_Z   = round(buttonpnt(1,2)); if slice_Z==0, slice_Z=1; end
                  %             slice_Y   = round(buttonpnt(1,1)); if slice_Y==0, slice_Y=1; end
                  %
                  %             % display new slice & update slice/slider values
                  %             if (1 <= slice_Y && slice_Y <= max_Y) && ...
                  %                     (1 <= slice_Z && slice_Z <= max_Z)
                  %                 slice_X = round(str2num(get(app.slice_X, 'String')));
                  %                 % adjust for flipped images
                  %                 slice_Z = max_Z+1-slice_Z;
                  %                 str = sprintf('%3.3g', stat(slice_X, slice_Y, slice_Z));
                  %                 set(app.voxelValue, 'String', str);
                  %                 updateImages(handles, stat, '', slice_Y, slice_Z);
            
               otherwise
                  error('Unknown axes');
            end % end switch axes
        end
        
        function removeToolBarButtons(app)
            removeItems = ([{'Save Figure'}, {'New Figure'}, {'Open File'}, {'Print Figure'}, {'Link Plot'}, {'Open Property Inspector'}, {'Insert Colorbar'}]);
            for i = 1:size(removeItems,2)
               listOfElements = findall(gcf);
               element = findall(listOfElements,'ToolTipString',string(removeItems(i)));
               set(element,'Visible','off');
            end
        end
        
        function varargout = removeVoltage(app, varargin)
         try
            varargout = app.f.removeVoltage(app, varargin);
         catch E
            if strcmp('MATLAB:unassignedOutputs', E.identifier)
               varargout = {};
            else
               rethrow(E);
            end
        end
        end
                
        function run_tool(app)
            app.f.run_tool(app);
        end
        
        function scroll_axes(app, hObject, eventdata, handles)
            tseries = app.f.getCurrentVoltage(app);
            handles = updateSETFigure(handles, tseries);
            guidata(hObject, handles);
        end
        
        function time_slider_updated(app, hObject)
            app.f.time_slider_updated(app, hObject)
        end
        
        function updateTimeSlider(app)
            app.f.updateTimeSlider(app);
        end        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function SpikeExtractionTool_OpeningFcn(app)
            % prepare for user confirmation when GUI is closed
            % app.figure1.CloseRequestFcn = 'closeGUI'; % It is commented
            % because it does not work for App. It is supposed to, but
            % there's an error I haven't figured out.
            try
               close(1); % Running the app opens a new figure sometimes. Close it.
            catch
               nan;
            end
            
            % Instantiate frontEndFunctions to access the shared functions
            % between GUIDE and App
            app.f = frontEndFunctions('app');
            %% Generate user data object to store
            app.data.num_tseries = 0;
            app.data.tseries     = [];   % initialise cell array of time series
            app.data.curr_tseries= [];   % current time series
            app.data.tseries_str = [];   % initialise names of time series
            app.data.used_names  = [];   % list of names already used by user
            app.data.last_tseries= [];   % last time series - restore if current ts deleted
            app.data.tlim        = [];   % time limits
            app.data.vlim        = [];   % voltage limits
            app.data.curr_tool   = [];   % record of user's current tool
            app.data.last_tool   = [];   % record of user's last tool
            app.data.params      = getDefaultToolParams; % default params for all tools & implementation methods
            % data.guihandles  = hObject;              % copy of all SET gui handles
            app.data.zoomPercentage = [0 0]; % Records currently chosen zoom value
            app.data.displacementPercentage = [0 0.5]; % Records currently chosen displacement value
            
            app.data.last_dir    = pwd;  % where they opened gui from
            
            handles = toggleSETGUIstate(app,'off'); % switch everything off until data's loaded            
                        
            % from property inspector in guide, under WindowScrollWheelFcn
            % @(hObject,eventdata)SpikeExtractionTool('figure1_WindowScrollWheelFcn',hObject,eventdata,guidata(hObject))
            
            % Initialize options for debug and loading window
            app.options.loadingWindowOn = true;
            app.options.debugOption = 'semi';
            app.options.rescaleOption = 'at_end';
            
            % Initialize tooltips
            handles = setTooltips(handles);
            
            % UIWAIT makes SpikeExtractionTool wait for user response (see UIRESUME)
            % uiwait(app.figure1);
            
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

        % Menu selected function: aboutMenuItem
        function aboutMenuItem_Callback(app, event)
            app.f.aboutMenuItem();
        end

        % Menu selected function: accessVoltageMenuBar
        function accessVoltageMenuBar_Callback(app, event)
            % hObject    handle to accessVoltageMenuBar (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            access_voltage_Callback(app);
        end

        % Button pushed function: access_voltage
        function access_voltage_Callback(app, event)
            [tseries, ~, ~, ts_name] = app.f.getCurrentVoltage(app);
            assignin('base', title2Str(ts_name,1,1,'_'), tseries);
        end

        % Button pushed function: add_voltage_button
        function add_voltage_button_Callback(app, event)
            mouseWaitingFunction(app.figure1,@add_voltage,app);
        end

        % Value changed function: automatic_params
        function automatic_params_Callback(app, event)
            app.f.automatic_params(app, event.Source);
        end

        % Menu selected function: clearDifferentMenuBar
        function clearDifferentMenu_Callback(app, event)
            % hObject    handle to clearDifferentMenu (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
               if app.data.num_tseries == 1
                  displayErrorMsg('There is only one voltage loaded.');
                  return
               else
                  % Find position of current mouse
                  get(app.figure1,'Parent');
                  screenHandle = get(app.figure1,'Parent');
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
                     @(hObject,eventdata)SpikeExtractionTool('clearVoltages',hObject,event,mainData,vtc));
               end
        end

        % Button pushed function: clear_voltage_button
        function clear_voltage_button_Callback(app, event)
            [tseries, ts_num, data_type, ts_name] = app.f.getCurrentVoltage(app);
            response = userConfirmation(['Delete ' ts_name '?'],...
               'Clear current time series?');
            if strcmp(response,'No')
               return
            end
            
            mouseWaitingFunction(app.figure1,@removeVoltage,app); % Instead of removeVoltage(handles);
        end

        % Value changed function: curr_signal
        function curr_signal_Callback(app, event)
            app.f.curr_signal(event.Source, event, app);
        end

        % Menu selected function: exitMenuBar
        function exitMenuBar_Callback(app, event)
            % hObject    handle to exitMenuBar (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            closeGUI(hObject,eventdata);
        end

        % Window scroll wheel function: figure1
        function figure1_WindowScrollWheelFcn(app, event)
            % hObject    handle to figure1 (see GCBO)
            % eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
            %	VerticalScrollCount: signed integer indicating direction and number of clicks
            %	VerticalScrollAmount: number of lines scrolled for each click
            % handles    structure with handles and user data (see GUIDATA)
            
            % Mouse pointer back to normal.
            set(app.figure1, 'pointer', 'arrow')
        end

        % Menu selected function: helpMenuItem
        function helpMenuItem_Callback(app, event)
            % hObject    handle to helpMenuItem (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
               open('resources\SEThelp.pdf');
        end

        % Button pushed function: load_voltage_button
        function load_voltage_button_Callback(app, event)
            % Warn if there is user data
            okToClear = false;
            if app.data.num_tseries == 0
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
            
            mouseWaitingFunction(app.figure1,@load_voltage,app);
        end

        % Value changed function: method_list
        function method_list_Callback(app, event)
            % update the tooltips
            tooltip_name = {['method_list_' lower(app.method_list.Value)]};
            setTooltips(app, {'method_list'}, getTooltips(tooltip_name));
        end

        % Button pushed function: new_figure
        function new_figure_Callback(app, event)
            app.f.new_figure(app);
        end

        % Menu selected function: plotNewFigMenuBar
        function plotNewFigMenuBar_Callback(app, event)
            % hObject    handle to plotNewFigMenuBar (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            new_figure_Callback(hObject,eventdata,handles);
        end

        % Button pushed function: run_tool_button
        function run_tool_button_Callback(app, event)
            mouseWaitingFunction(app.figure1,@run_tool,app);
        end

        % Menu selected function: saveVoltageMenuBar
        function saveVoltageMenuBar_Callback(app, event)
            % hObject    handle to saveVoltageMenuBar (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            save_voltage_Callback(hObject, eventdata, handles);
        end

        % Button pushed function: save_voltage
        function save_voltage_Callback(app, event)
            app.f.save_voltage(app);
        end

        % Value changed function: scroll_axes_slider
        function scroll_axes_slider_Callback(app, event)
            mouseWaitingFunction(app.figure1, @scroll_axes, hObject, eventdata, handles);
        end

        % Button pushed function: set_tool_params
        function set_tool_params_Callback(app, event)
            app.f.set_tool_params(app);
        end

        % Value changed function: time_max
        function time_max_Callback(app, event)
           app.f.time_max(app, event.Source); 
        end

        % Value changed function: time_min
        function time_min_Callback(app, event)
           app.f.time_min(app);
        end

        % Value changed function: time_slider
        function time_slider_Callback(app, event)
            if ~app.f.haveUserData(app)
               return;
            end
            mouseWaitingFunction(app.figure1,@time_slider_updated,app, event.Source);
        end

        % Button pushed function: toggleZoomButton
        function toggleZoomButton_Callback(app, event)
            app.f.toggleZoomButton(app, event.Source);
        end

        % Value changed function: tool_list
        function tool_list_Callback(app, event)
            app.f.tool_list(event.Source, app);
        end

        % Value changed function: voltage_max
        function voltage_max_Callback(app, event)
            app.f.voltage_max(app, event.Source);
        end

        % Value changed function: voltage_min
        function voltage_min_Callback(app, event)
            app.f.voltage_min(app, event.Source);
        end

        % Value changed function: voltage_slider
        function voltage_slider_Callback(app, event)
            if ~app.f.haveUserData(app)
               return;
            end
            
            mouseWaitingFunction(app.figure1,@voltage_slider_updated, app, event.Source);
        end
        
        function voltage_slider_updated(app, hObject)
            app.f.voltage_slider_updated(app, hObject);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create figure1 and hide until all components are created
            app.figure1 = uifigure('Visible', 'off');
            app.figure1.Color = [1 1 1];
            app.figure1.Position = [942 257 802 637];
            app.figure1.Name = 'Spike Extraction Tool';
            app.figure1.WindowScrollWheelFcn = createCallbackFcn(app, @figure1_WindowScrollWheelFcn, true);

            % Create fileMenu
            app.fileMenu = uimenu(app.figure1);
            app.fileMenu.Text = 'File';

            % Create loadMenuBar
            app.loadMenuBar = uimenu(app.fileMenu);
            app.loadMenuBar.Text = 'Load Voltage';

            % Create addVoltageMenuBar
            app.addVoltageMenuBar = uimenu(app.fileMenu);
            app.addVoltageMenuBar.Enable = 'off';
            app.addVoltageMenuBar.Text = 'Add Voltage';

            % Create clearMenuBar
            app.clearMenuBar = uimenu(app.fileMenu);
            app.clearMenuBar.Enable = 'off';
            app.clearMenuBar.Text = 'Clear Voltage';

            % Create clearSelectedMenuBar
            app.clearSelectedMenuBar = uimenu(app.clearMenuBar);
            app.clearSelectedMenuBar.Text = 'Clear Selected Voltage';

            % Create clearDifferentMenuBar
            app.clearDifferentMenuBar = uimenu(app.clearMenuBar);
            app.clearDifferentMenuBar.MenuSelectedFcn = createCallbackFcn(app, @clearDifferentMenu_Callback, true);
            app.clearDifferentMenuBar.Text = 'Clear a Different Voltage';

            % Create plotNewFigMenuBar
            app.plotNewFigMenuBar = uimenu(app.fileMenu);
            app.plotNewFigMenuBar.MenuSelectedFcn = createCallbackFcn(app, @plotNewFigMenuBar_Callback, true);
            app.plotNewFigMenuBar.Separator = 'on';
            app.plotNewFigMenuBar.Accelerator = 'N';
            app.plotNewFigMenuBar.Text = 'Plot in New Figure';

            % Create accessVoltageMenuBar
            app.accessVoltageMenuBar = uimenu(app.fileMenu);
            app.accessVoltageMenuBar.MenuSelectedFcn = createCallbackFcn(app, @accessVoltageMenuBar_Callback, true);
            app.accessVoltageMenuBar.Enable = 'off';
            app.accessVoltageMenuBar.Text = 'Access Voltage';

            % Create saveVoltageMenuBar
            app.saveVoltageMenuBar = uimenu(app.fileMenu);
            app.saveVoltageMenuBar.MenuSelectedFcn = createCallbackFcn(app, @saveVoltageMenuBar_Callback, true);
            app.saveVoltageMenuBar.Enable = 'off';
            app.saveVoltageMenuBar.Accelerator = 'S';
            app.saveVoltageMenuBar.Text = 'Save Voltage';

            % Create exitMenuBar
            app.exitMenuBar = uimenu(app.fileMenu);
            app.exitMenuBar.MenuSelectedFcn = createCallbackFcn(app, @exitMenuBar_Callback, true);
            app.exitMenuBar.Separator = 'on';
            app.exitMenuBar.Text = 'Exit';

            % Create helpMenu
            app.helpMenu = uimenu(app.figure1);
            app.helpMenu.Text = 'Help';

            % Create helpMenuItem
            app.helpMenuItem = uimenu(app.helpMenu);
            app.helpMenuItem.MenuSelectedFcn = createCallbackFcn(app, @helpMenuItem_Callback, true);
            app.helpMenuItem.Text = 'Help';

            % Create aboutMenuItem
            app.aboutMenuItem = uimenu(app.helpMenu);
            app.aboutMenuItem.MenuSelectedFcn = createCallbackFcn(app, @aboutMenuItem_Callback, true);
            app.aboutMenuItem.Separator = 'on';
            app.aboutMenuItem.Text = 'About';

            % Create load_panel
            app.load_panel = uipanel(app.figure1);
            app.load_panel.BackgroundColor = [1 0.9 1];
            app.load_panel.FontSize = 13;
            app.load_panel.Position = [20 472 193 125];

            % Create load_voltage_button
            app.load_voltage_button = uibutton(app.load_panel, 'push');
            app.load_voltage_button.ButtonPushedFcn = createCallbackFcn(app, @load_voltage_button_Callback, true);
            app.load_voltage_button.BusyAction = 'cancel';
            app.load_voltage_button.BackgroundColor = [0.9 0.7 1];
            app.load_voltage_button.FontSize = 13;
            app.load_voltage_button.FontWeight = 'bold';
            app.load_voltage_button.Position = [37 57 116 25];
            app.load_voltage_button.Text = 'Load voltage';

            % Create select_dataset
            app.select_dataset = uilabel(app.load_panel);
            app.select_dataset.BackgroundColor = [1 0.9 1];
            app.select_dataset.HorizontalAlignment = 'center';
            app.select_dataset.VerticalAlignment = 'top';
            app.select_dataset.FontSize = 19;
            app.select_dataset.FontWeight = 'bold';
            app.select_dataset.Position = [33 91 124 25];
            app.select_dataset.Text = 'Select dataset';

            % Create add_voltage_button
            app.add_voltage_button = uibutton(app.load_panel, 'push');
            app.add_voltage_button.ButtonPushedFcn = createCallbackFcn(app, @add_voltage_button_Callback, true);
            app.add_voltage_button.BusyAction = 'cancel';
            app.add_voltage_button.BackgroundColor = [0.9 0.7 1];
            app.add_voltage_button.FontSize = 13;
            app.add_voltage_button.FontWeight = 'bold';
            app.add_voltage_button.Position = [37 20 116 25];
            app.add_voltage_button.Text = 'Add voltage';

            % Create voltage_panel
            app.voltage_panel = uipanel(app.figure1);
            app.voltage_panel.BackgroundColor = [0.8 0.9 1];
            app.voltage_panel.FontSize = 13;
            app.voltage_panel.Position = [20 322 193 129];

            % Create display_voltage
            app.display_voltage = uilabel(app.voltage_panel);
            app.display_voltage.BackgroundColor = [0.8 0.9 1];
            app.display_voltage.HorizontalAlignment = 'center';
            app.display_voltage.VerticalAlignment = 'top';
            app.display_voltage.FontSize = 19;
            app.display_voltage.FontWeight = 'bold';
            app.display_voltage.Position = [34 92 124 26];
            app.display_voltage.Text = 'Display voltage';

            % Create curr_signal
            app.curr_signal = uidropdown(app.voltage_panel);
            app.curr_signal.Items = {};
            app.curr_signal.ValueChangedFcn = createCallbackFcn(app, @curr_signal_Callback, true);
            app.curr_signal.BusyAction = 'cancel';
            app.curr_signal.Interruptible = 'off';
            app.curr_signal.FontSize = 13;
            app.curr_signal.FontWeight = 'bold';
            app.curr_signal.Position = [33 54 125 30];
            app.curr_signal.Value = {};

            % Create clear_voltage_button
            app.clear_voltage_button = uibutton(app.voltage_panel, 'push');
            app.clear_voltage_button.ButtonPushedFcn = createCallbackFcn(app, @clear_voltage_button_Callback, true);
            app.clear_voltage_button.BusyAction = 'cancel';
            app.clear_voltage_button.BackgroundColor = [0.2 0.6 1];
            app.clear_voltage_button.FontSize = 13;
            app.clear_voltage_button.FontWeight = 'bold';
            app.clear_voltage_button.Position = [38 21 116 25];
            app.clear_voltage_button.Text = 'Clear voltage';

            % Create tool_panel
            app.tool_panel = uipanel(app.figure1);
            app.tool_panel.BackgroundColor = [1 0.9 0.6];
            app.tool_panel.FontSize = 13;
            app.tool_panel.Position = [524 8 272 198];

            % Create tools_title
            app.tools_title = uilabel(app.tool_panel);
            app.tools_title.BackgroundColor = [1 0.9 0.6];
            app.tools_title.HorizontalAlignment = 'center';
            app.tools_title.VerticalAlignment = 'top';
            app.tools_title.FontSize = 21;
            app.tools_title.FontWeight = 'bold';
            app.tools_title.Position = [29 160 220 31];
            app.tools_title.Text = 'Available tools';

            % Create tool_list
            app.tool_list = uidropdown(app.tool_panel);
            app.tool_list.Items = {'Rescale', 'Denoise', 'Identify AP templates', 'Extract spikes', 'Utilities'};
            app.tool_list.ValueChangedFcn = createCallbackFcn(app, @tool_list_Callback, true);
            app.tool_list.FontSize = 13;
            app.tool_list.FontWeight = 'bold';
            app.tool_list.BackgroundColor = [1 1 1];
            app.tool_list.Position = [125 128 141 27];
            app.tool_list.Value = 'Rescale';

            % Create select_tool
            app.select_tool = uilabel(app.tool_panel);
            app.select_tool.BackgroundColor = [0.93 0.69 0.13];
            app.select_tool.HorizontalAlignment = 'center';
            app.select_tool.VerticalAlignment = 'top';
            app.select_tool.FontSize = 13;
            app.select_tool.FontWeight = 'bold';
            app.select_tool.Position = [9 133 103 18];
            app.select_tool.Text = 'Select tool';

            % Create method_list
            app.method_list = uidropdown(app.tool_panel);
            app.method_list.Items = {'Particle filter', 'Recursive least squares', 'Recursive mean', 'Variance'};
            app.method_list.ValueChangedFcn = createCallbackFcn(app, @method_list_Callback, true);
            app.method_list.FontSize = 13;
            app.method_list.FontWeight = 'bold';
            app.method_list.BackgroundColor = [1 1 1];
            app.method_list.Position = [125 86 141 27];
            app.method_list.Value = 'Particle filter';

            % Create select_method
            app.select_method = uilabel(app.tool_panel);
            app.select_method.BackgroundColor = [0.93 0.69 0.13];
            app.select_method.HorizontalAlignment = 'center';
            app.select_method.VerticalAlignment = 'top';
            app.select_method.FontSize = 13;
            app.select_method.FontWeight = 'bold';
            app.select_method.Position = [9 90 105 19];
            app.select_method.Text = 'Select method';

            % Create set_tool_params
            app.set_tool_params = uibutton(app.tool_panel, 'push');
            app.set_tool_params.ButtonPushedFcn = createCallbackFcn(app, @set_tool_params_Callback, true);
            app.set_tool_params.BusyAction = 'cancel';
            app.set_tool_params.BackgroundColor = [0.93 0.69 0.13];
            app.set_tool_params.FontSize = 13;
            app.set_tool_params.FontWeight = 'bold';
            app.set_tool_params.Position = [9 38 106 30];
            app.set_tool_params.Text = 'Set parameters';

            % Create run_tool_button
            app.run_tool_button = uibutton(app.tool_panel, 'push');
            app.run_tool_button.ButtonPushedFcn = createCallbackFcn(app, @run_tool_button_Callback, true);
            app.run_tool_button.BusyAction = 'cancel';
            app.run_tool_button.Interruptible = 'off';
            app.run_tool_button.BackgroundColor = [0.93 0.69 0.13];
            app.run_tool_button.FontSize = 13;
            app.run_tool_button.FontWeight = 'bold';
            app.run_tool_button.Position = [158 39 105 30];
            app.run_tool_button.Text = 'Run tool';

            % Create automatic_params
            app.automatic_params = uicheckbox(app.tool_panel);
            app.automatic_params.ValueChangedFcn = createCallbackFcn(app, @automatic_params_Callback, true);
            app.automatic_params.Tooltip = 'sprintf(''If ticked, the parameters are set autmatically. \nUntick this box to set parameters manually.'')';
            app.automatic_params.Text = 'Automatic parameters';
            app.automatic_params.FontSize = 13;
            app.automatic_params.FontWeight = 'bold';
            app.automatic_params.Position = [10 8 173 23];

            % Create param_panel
            app.param_panel = uipanel(app.figure1);
            app.param_panel.BackgroundColor = [0.8 1 0.6];
            app.param_panel.FontSize = 13;
            app.param_panel.Position = [39 9 448 198];

            % Create gen_params
            app.gen_params = uilabel(app.param_panel);
            app.gen_params.BackgroundColor = [0.8 1 0.6];
            app.gen_params.VerticalAlignment = 'top';
            app.gen_params.FontSize = 21;
            app.gen_params.FontWeight = 'bold';
            app.gen_params.Position = [9 163 431 27];
            app.gen_params.Text = 'Generation parameters';

            % Create param1
            app.param1 = uilabel(app.param_panel);
            app.param1.BackgroundColor = [0.2 0.67 0.3];
            app.param1.HorizontalAlignment = 'center';
            app.param1.VerticalAlignment = 'top';
            app.param1.FontSize = 16;
            app.param1.FontWeight = 'bold';
            app.param1.Visible = 'off';
            app.param1.Position = [8 128 99 27];
            app.param1.Text = 'off';

            % Create val1
            app.val1 = uilabel(app.param_panel);
            app.val1.BackgroundColor = [0.2 0.67 0.3];
            app.val1.HorizontalAlignment = 'center';
            app.val1.VerticalAlignment = 'top';
            app.val1.FontSize = 16;
            app.val1.FontWeight = 'bold';
            app.val1.Visible = 'off';
            app.val1.Position = [117 128 99 27];
            app.val1.Text = 'off';

            % Create param2
            app.param2 = uilabel(app.param_panel);
            app.param2.BackgroundColor = [0.2 0.67 0.3];
            app.param2.HorizontalAlignment = 'center';
            app.param2.VerticalAlignment = 'top';
            app.param2.FontSize = 16;
            app.param2.FontWeight = 'bold';
            app.param2.Visible = 'off';
            app.param2.Position = [8 89 99 27];
            app.param2.Text = 'off';

            % Create val2
            app.val2 = uilabel(app.param_panel);
            app.val2.BackgroundColor = [0.2 0.67 0.3];
            app.val2.HorizontalAlignment = 'center';
            app.val2.VerticalAlignment = 'top';
            app.val2.FontSize = 16;
            app.val2.FontWeight = 'bold';
            app.val2.Visible = 'off';
            app.val2.Position = [117 89 99 27];
            app.val2.Text = 'off';

            % Create param3
            app.param3 = uilabel(app.param_panel);
            app.param3.BackgroundColor = [0.2 0.67 0.3];
            app.param3.HorizontalAlignment = 'center';
            app.param3.VerticalAlignment = 'top';
            app.param3.FontSize = 16;
            app.param3.FontWeight = 'bold';
            app.param3.Visible = 'off';
            app.param3.Position = [8 50 99 27];
            app.param3.Text = 'off';

            % Create val3
            app.val3 = uilabel(app.param_panel);
            app.val3.BackgroundColor = [0.2 0.67 0.3];
            app.val3.HorizontalAlignment = 'center';
            app.val3.VerticalAlignment = 'top';
            app.val3.FontSize = 16;
            app.val3.FontWeight = 'bold';
            app.val3.Visible = 'off';
            app.val3.Position = [117 50 99 27];
            app.val3.Text = 'off';

            % Create param4
            app.param4 = uilabel(app.param_panel);
            app.param4.BackgroundColor = [0.2 0.67 0.3];
            app.param4.HorizontalAlignment = 'center';
            app.param4.VerticalAlignment = 'top';
            app.param4.FontSize = 16;
            app.param4.FontWeight = 'bold';
            app.param4.Visible = 'off';
            app.param4.Position = [8 11 99 27];
            app.param4.Text = 'off';

            % Create val4
            app.val4 = uilabel(app.param_panel);
            app.val4.BackgroundColor = [0.2 0.67 0.3];
            app.val4.HorizontalAlignment = 'center';
            app.val4.VerticalAlignment = 'top';
            app.val4.FontSize = 16;
            app.val4.FontWeight = 'bold';
            app.val4.Visible = 'off';
            app.val4.Position = [117 11 99 27];
            app.val4.Text = 'off';

            % Create param5
            app.param5 = uilabel(app.param_panel);
            app.param5.BackgroundColor = [0.2 0.67 0.3];
            app.param5.HorizontalAlignment = 'center';
            app.param5.VerticalAlignment = 'top';
            app.param5.FontSize = 16;
            app.param5.FontWeight = 'bold';
            app.param5.Visible = 'off';
            app.param5.Position = [232 128 99 27];
            app.param5.Text = 'off';

            % Create val5
            app.val5 = uilabel(app.param_panel);
            app.val5.BackgroundColor = [0.2 0.67 0.3];
            app.val5.HorizontalAlignment = 'center';
            app.val5.VerticalAlignment = 'top';
            app.val5.FontSize = 16;
            app.val5.FontWeight = 'bold';
            app.val5.Visible = 'off';
            app.val5.Position = [341 128 99 27];
            app.val5.Text = 'off';

            % Create param6
            app.param6 = uilabel(app.param_panel);
            app.param6.BackgroundColor = [0.2 0.67 0.3];
            app.param6.HorizontalAlignment = 'center';
            app.param6.VerticalAlignment = 'top';
            app.param6.FontSize = 16;
            app.param6.FontWeight = 'bold';
            app.param6.Visible = 'off';
            app.param6.Position = [232 89 99 27];
            app.param6.Text = 'off';

            % Create val6
            app.val6 = uilabel(app.param_panel);
            app.val6.BackgroundColor = [0.2 0.67 0.3];
            app.val6.HorizontalAlignment = 'center';
            app.val6.VerticalAlignment = 'top';
            app.val6.FontSize = 16;
            app.val6.FontWeight = 'bold';
            app.val6.Visible = 'off';
            app.val6.Position = [340 89 99 27];
            app.val6.Text = 'off';

            % Create param7
            app.param7 = uilabel(app.param_panel);
            app.param7.BackgroundColor = [0.2 0.67 0.3];
            app.param7.HorizontalAlignment = 'center';
            app.param7.VerticalAlignment = 'top';
            app.param7.FontSize = 16;
            app.param7.FontWeight = 'bold';
            app.param7.Visible = 'off';
            app.param7.Position = [232 50 99 27];
            app.param7.Text = 'off';

            % Create val7
            app.val7 = uilabel(app.param_panel);
            app.val7.BackgroundColor = [0.2 0.67 0.3];
            app.val7.HorizontalAlignment = 'center';
            app.val7.VerticalAlignment = 'top';
            app.val7.FontSize = 16;
            app.val7.FontWeight = 'bold';
            app.val7.Visible = 'off';
            app.val7.Position = [340 50 99 27];
            app.val7.Text = 'off';

            % Create param8
            app.param8 = uilabel(app.param_panel);
            app.param8.BackgroundColor = [0.2 0.67 0.3];
            app.param8.HorizontalAlignment = 'center';
            app.param8.VerticalAlignment = 'top';
            app.param8.FontSize = 16;
            app.param8.FontWeight = 'bold';
            app.param8.Visible = 'off';
            app.param8.Position = [232 11 99 27];
            app.param8.Text = 'off';

            % Create val8
            app.val8 = uilabel(app.param_panel);
            app.val8.BackgroundColor = [0.2 0.67 0.3];
            app.val8.HorizontalAlignment = 'center';
            app.val8.VerticalAlignment = 'top';
            app.val8.FontSize = 16;
            app.val8.FontWeight = 'bold';
            app.val8.Visible = 'off';
            app.val8.Position = [340 11 99 27];
            app.val8.Text = 'off';

            % Create save_voltage
            app.save_voltage = uibutton(app.figure1, 'push');
            app.save_voltage.ButtonPushedFcn = createCallbackFcn(app, @save_voltage_Callback, true);
            app.save_voltage.BackgroundColor = [1 0.6 0.6];
            app.save_voltage.FontSize = 13;
            app.save_voltage.FontWeight = 'bold';
            app.save_voltage.Position = [603 605 133 26];
            app.save_voltage.Text = 'Save voltage';

            % Create plot_panel
            app.plot_panel = uipanel(app.figure1);
            app.plot_panel.BorderType = 'none';
            app.plot_panel.BackgroundColor = [1 1 1];
            app.plot_panel.FontWeight = 'bold';
            app.plot_panel.FontSize = 21;
            app.plot_panel.Position = [234 261 545 334];

            % Create axes1
            app.axes1 = uiaxes(app.plot_panel);
            app.axes1.FontSize = 13;
            app.axes1.XLim = [0.000199999994947575 600.000183105469];
            app.axes1.YLim = [-5 4.578857421875];
            app.axes1.ColorOrderIndex = 2;
            app.axes1.Box = 'on';
            app.axes1.NextPlot = 'replace';
            app.axes1.BackgroundColor = [1 1 1];
            app.axes1.Position = [54 10 452 300];

            % Create voltage_slider
            app.voltage_slider = uislider(app.plot_panel);
            app.voltage_slider.Limits = [0 1];
            app.voltage_slider.MajorTicks = [];
            app.voltage_slider.Orientation = 'vertical';
            app.voltage_slider.ValueChangedFcn = createCallbackFcn(app, @voltage_slider_Callback, true);
            app.voltage_slider.MinorTicks = [];
            app.voltage_slider.BusyAction = 'cancel';
            app.voltage_slider.FontSize = 13;
            app.voltage_slider.FontColor = [0.64 0.08 0.18];
            app.voltage_slider.Position = [16 91 3 188];

            % Create voltage_max
            app.voltage_max = uieditfield(app.plot_panel, 'text');
            app.voltage_max.ValueChangedFcn = createCallbackFcn(app, @voltage_max_Callback, true);
            app.voltage_max.BusyAction = 'cancel';
            app.voltage_max.HorizontalAlignment = 'center';
            app.voltage_max.FontSize = 13;
            app.voltage_max.BackgroundColor = [1 0.6 0.6];
            app.voltage_max.Position = [0 289 44 23];
            app.voltage_max.Value = '4.5789';

            % Create voltage_min
            app.voltage_min = uieditfield(app.plot_panel, 'text');
            app.voltage_min.ValueChangedFcn = createCallbackFcn(app, @voltage_min_Callback, true);
            app.voltage_min.BusyAction = 'cancel';
            app.voltage_min.HorizontalAlignment = 'center';
            app.voltage_min.FontSize = 13;
            app.voltage_min.BackgroundColor = [1 0.6 0.6];
            app.voltage_min.Position = [0 55 44 23];
            app.voltage_min.Value = '-5';

            % Create time_slider
            app.time_slider = uislider(app.figure1);
            app.time_slider.Limits = [0 1];
            app.time_slider.MajorTicks = [];
            app.time_slider.ValueChangedFcn = createCallbackFcn(app, @time_slider_Callback, true);
            app.time_slider.MinorTicks = [];
            app.time_slider.FontSize = 13;
            app.time_slider.Position = [405 241 251 3];

            % Create time_max
            app.time_max = uieditfield(app.figure1, 'text');
            app.time_max.ValueChangedFcn = createCallbackFcn(app, @time_max_Callback, true);
            app.time_max.HorizontalAlignment = 'center';
            app.time_max.FontSize = 13;
            app.time_max.BackgroundColor = [1 0.6 0.6];
            app.time_max.Position = [662 237 64 25];
            app.time_max.Value = '600.0002';

            % Create time_min
            app.time_min = uieditfield(app.figure1, 'text');
            app.time_min.ValueChangedFcn = createCallbackFcn(app, @time_min_Callback, true);
            app.time_min.HorizontalAlignment = 'center';
            app.time_min.FontSize = 13;
            app.time_min.BackgroundColor = [1 0.6 0.6];
            app.time_min.Position = [333 237 64 25];
            app.time_min.Value = '0.0002';

            % Create access_voltage
            app.access_voltage = uibutton(app.figure1, 'push');
            app.access_voltage.ButtonPushedFcn = createCallbackFcn(app, @access_voltage_Callback, true);
            app.access_voltage.BackgroundColor = [1 0.6 0.6];
            app.access_voltage.FontSize = 13;
            app.access_voltage.FontWeight = 'bold';
            app.access_voltage.Position = [454 605 133 26];
            app.access_voltage.Text = 'Access voltage';

            % Create new_figure
            app.new_figure = uibutton(app.figure1, 'push');
            app.new_figure.ButtonPushedFcn = createCallbackFcn(app, @new_figure_Callback, true);
            app.new_figure.BackgroundColor = [1 0.6 0.6];
            app.new_figure.FontSize = 13;
            app.new_figure.FontWeight = 'bold';
            app.new_figure.Position = [306 605 133 26];
            app.new_figure.Text = 'Plot in new figure';

            % Create scroll_axes_slider
            app.scroll_axes_slider = uislider(app.figure1);
            app.scroll_axes_slider.Limits = [0 1];
            app.scroll_axes_slider.MajorTicks = [];
            app.scroll_axes_slider.Orientation = 'vertical';
            app.scroll_axes_slider.ValueChangedFcn = createCallbackFcn(app, @scroll_axes_slider_Callback, true);
            app.scroll_axes_slider.MinorTicks = [];
            app.scroll_axes_slider.Visible = 'off';
            app.scroll_axes_slider.FontSize = 13;
            app.scroll_axes_slider.FontColor = [0.64 0.08 0.18];
            app.scroll_axes_slider.Position = [781 318 3 250];
            app.scroll_axes_slider.Value = 1;

            % Create toggleZoomButton
            app.toggleZoomButton = uibutton(app.figure1, 'push');
            app.toggleZoomButton.ButtonPushedFcn = createCallbackFcn(app, @toggleZoomButton_Callback, true);
            app.toggleZoomButton.Icon = 'fig/toggleZoomButton_image.png';
            app.toggleZoomButton.IconAlignment = 'center';
            app.toggleZoomButton.BackgroundColor = [1 0.6 0.6];
            app.toggleZoomButton.FontSize = 11;
            app.toggleZoomButton.Position = [286 233 32 32];
            app.toggleZoomButton.Text = '';
            app.toggleZoomButton.UserData = 'zoom';

            % Show the figure after all components are created
            app.figure1.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SpikeExtractionApp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.figure1)

            % Execute the startup function
            runStartupFcn(app, @SpikeExtractionTool_OpeningFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)            
            % Delete UIFigure when app is deleted
            delete(app.figure1)
        end
    end
end