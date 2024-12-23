function [method_params,cancel] = requestUserParamConfig_identify_AP_template(method_params,dlg_name)
    % requestUserParamConfig_identify_AP_template Description: This
    % function makes a special UI window for the identifyAPtemplate tool in
    % spikeExtractionTool
    % ______________
    %[method_params,cancel] = requestUserParamConfig_identify_AP_template(method_params,dlg_name)
    % Inputs:
    %   - method_params: a struct containing a field for each changable
    %   parameter. each field contains a struct with fields value, name
    %   descript,type, and units
    %   - dlg_name: name for the dialog box
    % Outputs:
    %   - method_params: value is changed by user for each field
    %   - cancel: if the user pressed cancel - return this
    % ______________
    % Author: Kathrine Clarke 
    % Date: 01-Dec-2023

    % Output help if there are no inputs
    if nargin == 0
       help requestUserParamConfig_identify_AP_template
       return;
    end

    % Initialize data structure to store parameters
    def_parameters = struct('N_phases', 2, 'V_first', 1, 'A_phase', 1, 'Vthresh', [3 3], 'Phase_width', [1.2 1.3]);
    parameters = repmat(def_parameters,0,1);
    objects = struct();
    
    % Create the main figure
    fig = figure('Name', 'Spike Settings', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 500]);

    % Panel 1: Radio buttons, 'Add', 'Delete'
    panel1 = uipanel('Parent', fig, 'Title', 'Spike Options', 'Position', [0.02, 0.1, 0.25, 0.8]);
    
    objects.radioGroup = uibuttongroup('Parent', panel1, 'Position', [0.1, 0.65, 0.8, 0.3], 'Title', 'Select Spike');
    objects.addBtn = uicontrol('Parent', panel1, 'Style', 'pushbutton', 'String', 'Add', 'Position', [10, 10, 60, 25], 'Callback', {@addSpike, parameters, objects,def_parameters});
    objects.deleteBtn = uicontrol('Parent', panel1, 'Style', 'pushbutton', 'String', 'Delete', 'Position', [80, 10, 60, 25], 'Callback', {@deleteSpike, parameters, objects});
    
    % Panel 2: Integer input boxes
    panel2 = uipanel('Parent', fig, 'Title', 'Phase Settings', 'Position', [0.3, 0.1, 0.25, 0.8]);

    objects.numPhasesLabel = uicontrol('Parent', panel2, 'Style', 'text', 'String', 'Number of Phases:', 'Position', [10, 90, 120, 20]);
    objects.voltageLabel = uicontrol('Parent', panel2, 'Style', 'text', 'String', 'Voltage of First Phase:', 'Position', [10, 60, 150, 20]);
    objects.alignmentLabel = uicontrol('Parent', panel2, 'Style', 'text', 'String', 'Alignment Phase:', 'Position', [10, 30, 120, 20]);

    objects.numPhasesEdit = uicontrol('Parent', panel2, 'Style', 'edit', 'Position', [160, 90, 50, 20], 'Callback', {@updatePanel3, parameters, objects});
    objects.voltageEdit = uicontrol('Parent', panel2, 'Style', 'edit', 'Position', [160, 60, 50, 20], 'Callback', {@updatePanel3, parameters, objects});
    objects.alignmentEdit = uicontrol('Parent', panel2, 'Style', 'edit', 'Position', [160, 30, 50, 20], 'Callback', {@updatePanel3, parameters, objects});

    % Panel 3: User input controls
    panel3 = uipanel('Parent', fig, 'Title', 'User Inputs', 'Position', [0.6, 0.1, 0.35, 0.8]);

    % Axes for plotting
    objects.ax = axes('Parent', fig, 'Position', [0.6, 0.1, 0.35, 0.8]);

    % 'Done' and 'Cancel' buttons
    objects.doneBtn = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Done', 'Position', [900, 20, 70, 30], 'Callback', @doneCallback);
    objects.cancelBtn = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Position', [980, 20, 70, 30], 'Callback', @cancelCallback);
    objects.addBtn = uicontrol('Parent', panel1, 'Style', 'pushbutton', 'String', 'Add', 'Position', [10, 10, 60, 25], 'Callback', {@addSpike, parameters, objects,def_parameters});
    
    [parameters, objects] = addSpike([],[],parameters, objects,def_parameters);
    
    % Functions
    function [parameters,objects] = addSpike(~, ~, parameters, objects, def_parameters)
        spikeNumber = numel(parameters) + 1;
        parameters(spikeNumber) = def_parameters;

        % Create a new radio button for the added spike
        objects.spike_buttons(spikeNumber) = uicontrol('Parent', objects.radioGroup, 'Style', 'radiobutton', 'String', ['Spike ', num2str(spikeNumber)], 'UserData', spikeNumber, 'Position', [10, 10 + 25 * (spikeNumber-1), 120, 20], 'Callback', {@updatePanels,parameters,objects});

        objects = updatePanels([],[],parameters,objects);
        objects = updatePanel3([], [], parameters, objects);
    end

    function deleteSpike(~, ~, parameters, objects)
        selectedIdx = get(get(objects.radioGroup, 'SelectedObject'), 'UserData');
        if ~isempty(selectedIdx)
            delete(objects.radioGroup.SelectedObject);
            parameters(selectedIdx) = [];
            objects = updatePanels();
            objects = updatePanel3([], [], parameters, objects);
        end
    end

    function objects = updatePanels(~, ~, parameters,objects)
        selectedIdx = get(get(objects.radioGroup, 'SelectedObject'), 'UserData');
        if ~isempty(selectedIdx)
            selectedSpike = parameters(selectedIdx);
            set(objects.numPhasesEdit, 'String', num2str(selectedSpike.N_phases));
            set(objects.voltageEdit, 'String', num2str(selectedSpike.V_first));
            set(objects.alignmentEdit, 'String', num2str(selectedSpike.A_phase));
        end
    end

    function objects = updatePanel3(~, ~, parameters, objects)
        selectedIdx = get(get(objects.radioGroup, 'SelectedObject'), 'UserData');
        if ~isempty(selectedIdx)
            numPhases = str2double(get(objects.numPhasesEdit, 'String'));
            voltage = str2double(get(objects.voltageEdit, 'String'));
            alignment = str2double(get(objects.alignmentEdit, 'String'));

            parameters(selectedIdx).N_phases = numPhases;
            parameters(selectedIdx).V_first = voltage;
            parameters(selectedIdx).A_phase = alignment;

            % Clear existing controls in Panel 3
            delete(findobj(panel3, 'Type', 'uicontrol'));

            % Create float input boxes dynamically in Panel 3
            for i = 1:(2 * numPhases)
                uicontrol('Parent', panel3, 'Style', 'edit', 'Position', [10, 90 - i * 30, 50, 20], 'String', '0.0', 'Callback', {@updatePhaseControls, parameters, selectedIdx, objects});
            end

            % Call the make_plot function
            make_plot(parameters, objects);
        end
    end

    function updatePhaseControls(~, ~, parameters, selectedIdx, objects)
        numPhases = parameters(selectedIdx).N_phases;
        Vthresh = zeros(numPhases, 1);
        Phase_width = zeros(numPhases, 1);

        for i = 1:(2 * numPhases)
            value = str2double(get(findobj(panel3, 'Style', 'edit', 'Position', [10, 90 - i * 30, 50, 20]), 'String'));
            if i <= numPhases
                Vthresh(i) = value;
            else
                Phase_width(i - numPhases) = value;
            end
        end

        parameters(selectedIdx).Vthresh = Vthresh;
        parameters(selectedIdx).Phase_width = Phase_width;

        % Call the make_plot function
        make_plot(parameters, objects);
    end

    function doneCallback(~, ~)
        close(fig);
    end

    function cancelCallback(~, ~)
        parameters = [];
        close(fig);
    end
end

% Function to make the plot based on user input controls
function make_plot(parameters, objects)
    % Extract parameters for plotting
    % You can customize this part based on your plotting requirements
    disp('Make your plot here using parameters:');
    disp(parameters);
end

