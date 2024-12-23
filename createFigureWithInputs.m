function parameters = createFigureWithInputs(defaultSpikes)
    % Check if defaultSpikes is provided, otherwise use default values
    if nargin < 1
        defaultSpikes = struct('N_phases', 2, 'V_first', 1, 'A_phase', 1, 'Vthresh', [3, 3], 'Phase_width', [1.3, 1.2]);
    end

    % Initialize data structure to store parameters
    parameters = defaultSpikes;

    % Create the main figure
    objects.fig = figure('Name', 'Spike Settings', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);

    % Panel 2: Integer input boxes
    panel2 = uipanel('Parent', objects.fig, 'Title', 'Phase Settings', 'Position', [0.1, 0.1, 0.4, 0.8]);

    objects.numPhasesLabel = uicontrol('Parent', panel2, 'Style', 'text', 'String', 'Number of Phases:', 'Position', [10, 90, 120, 20]);
    objects.voltageLabel = uicontrol('Parent', panel2, 'Style', 'text', 'String', 'Voltage of First Phase:', 'Position', [10, 60, 150, 20]);
    objects.alignmentLabel = uicontrol('Parent', panel2, 'Style', 'text', 'String', 'Alignment Phase:', 'Position', [10, 30, 120, 20]);

    objects.numPhasesEdit = uicontrol('Parent', panel2, 'Style', 'edit', 'String', num2str(defaultSpikes.N_phases), 'Position', [160, 90, 50, 20]);
    objects.voltageEdit = uicontrol('Parent', panel2, 'Style', 'edit', 'String', num2str(defaultSpikes.V_first), 'Position', [160, 60, 50, 20]);
    objects.alignmentEdit = uicontrol('Parent', panel2, 'Style', 'edit', 'String', num2str(defaultSpikes.A_phase), 'Position', [160, 30, 50, 20]);

    % Panel 3: User input controls
    panel3 = uipanel('Parent', objects.fig, 'Title', 'User Inputs', 'Position', [0.5, 0.1, 0.4, 0.8]);

    objects.vthreshLabel = uicontrol('Parent', panel3, 'Style', 'text', 'String', 'Voltage Thresholds (stds):', 'Position', [10, 90, 120, 20]);
    objects.phasetimeLabel = uicontrol('Parent', panel3, 'Style', 'text', 'String', 'Duration of phase (ms):', 'Position', [10, 60, 150, 20]);

    objects.vthreshEdit = uicontrol('Parent', panel3, 'Style', 'edit', 'String', num2str(defaultSpikes.Vthresh), 'Position', [160, 90, 100, 20]);
    objects.phasetimeEdit = uicontrol('Parent', panel3, 'Style', 'edit', 'String', num2str(defaultSpikes.Phase_width), 'Position', [160, 60, 100, 20]);

    % Axes for plotting
    objects.ax = axes('Parent', panel3, 'Position', [0.1, 0.5, 0.8, 0.4]);

    % 'Done' and 'Cancel' buttons
    objects.doneBtn = uicontrol('Parent', objects.fig, 'Style', 'pushbutton', 'String', 'Done', 'Position', [500, 20, 70, 30]);
    objects.cancelBtn = uicontrol('Parent', objects.fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Position', [580, 20, 70, 30]);

    % Set callbacks
    set(objects.numPhasesEdit, 'Callback', {@somethingChanged, parameters, objects});
    set(objects.voltageEdit, 'Callback', {@somethingChanged, parameters, objects});
    set(objects.alignmentEdit, 'Callback', {@somethingChanged, parameters, objects});
    set(objects.vthreshEdit, 'Callback', {@somethingChanged, parameters, objects});
    set(objects.phasetimeEdit, 'Callback', {@somethingChanged, parameters, objects});

    set(objects.doneBtn, 'Callback', @doneCallback);
    set(objects.cancelBtn, 'Callback', @cancelCallback);
    make_plot(parameters, objects)
    % Nested functions
    function somethingChanged(~, ~, ~, ~)
        updateParameters();
        make_plot(parameters, objects);
    end

    function updateParameters()
        Nphase_old = parameters.N_phases;
        parameters.N_phases = str2double(get(objects.numPhasesEdit, 'String'));
        parameters.V_first = str2double(get(objects.voltageEdit, 'String'));
        parameters.A_phase = str2double(get(objects.alignmentEdit, 'String'));
        parameters.Vthresh = str2num(get(objects.vthreshEdit, 'String'));
        parameters.Phase_width = str2num(get(objects.phasetimeEdit, 'String'));
        
        if Nphase_old >  parameters.N_phases
           parameters.Vthresh = parameters.Vthresh(1:parameters.N_phases);
           parameters.Phase_width = parameters.Phase_width(1:parameters.N_phases);
        end
        if Nphase_old <  parameters.N_phases
           parameters.Vthresh = [parameters.Vthresh 2*ones(1, parameters.N_phases-Nphase_old)];
           parameters.Phase_width = [parameters.Phase_width 2*ones(1, parameters.N_phases-Nphase_old)];
        end
        set(objects.vthreshEdit, 'String',num2str( parameters.Vthresh ));
        set(objects.phasetimeEdit, 'String',num2str(parameters.Phase_width));
        
    end

    function doneCallback(~, ~)
        close(objects.fig);
    end

    function cancelCallback(~, ~)
        close(objects.fig);
    end

    function make_plot(parameters, objects)
        axes(objects.ax);
        make_example_spike(parameters.N_phases, parameters.V_first, parameters.A_phase, parameters.Vthresh, parameters.Phase_width);
    end

    % Wait for the figure to close
    waitfor(objects.fig);
end
