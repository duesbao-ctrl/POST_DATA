classdef PostDataApp < handle
%POSTDATAAPP Full MATLAB R2016b-compatible POST_DATA desktop application.
% UI orchestration only; calculation, rendering, and export are services.

    properties
        Figure
        FileEdit
        TaskPopup
        TaskValues
        SelectionPopup
        SelectionEdit
        SlurmEdit
        SlurmBrowseButton
        ProgressPopup
        ResultLevelPopup
        RunButton
        CancelButton
        StatusText
        PlotPanel
        PlotAxes
        PlotViewPopup
        PlotViewIds
        PlotViewCountText
        OpenSeparateFiguresCheckbox
        PlotFigureHandles
        PlotFigureAxes
        PlotFigureViewIds
        TabGroup
        VisualizationTab
        PlotDataTab
        ComputeTab
        ComputeHeader
        ComputeTable
        ParameterNameText
        ParameterHelpText
        ParameterEdit
        ParameterPopup
        ParameterCheckbox
        ParameterApplyButton
        SelectedParameterIndex
        PlotTable
        PlotPresetPopup
        PlotParameterNameText
        PlotParameterHelpText
        PlotParameterEdit
        PlotParameterPopup
        PlotParameterCheckbox
        PlotParameterApplyButton
        SelectedPlotParameterIndex
        OutputTable
        ResultTable
        PlotDataViewPopup
        PlotDataViewIds
        PlotDataTable
        PlotDataInfoText
        PlotDataPageText
        PlotDataPreviousButton
        PlotDataNextButton
        PlotDataCopyHeadersCheckbox
        CurrentPlotData
        PlotDataPage
        PlotDataPageSize
        PlotDataSelectedRows
        CurrentResult
        OptionCatalog
        PlotCatalog
        OutputCatalog
        ActiveAnalysisType
        OptionDataCache
        CancelRequested
    end

    methods
        function obj = PostDataApp()
            obj.buildInterface();
            obj.resetAllOptions();
            obj.setStatus(pd_ui_text('Select data and analysis; validate parameters; run; review or export results'));
        end

        function buildInterface(obj)
            obj.Figure = figure('Name', pd_ui_text('POST_DATA2 MATLAB post-processing software'), ...
                'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'figure', ...
                'Tag', 'POST_DATA2_MainFigure', 'Color', [0.94, 0.94, 0.94], ...
                'Position', defaultWindowPosition());
            set(obj.Figure, 'CloseRequestFcn', @(src, evt)obj.closeApplication(src, evt));
            set(obj.Figure, 'WindowKeyPressFcn', @(src, evt)obj.handleShortcut(src, evt));
            obj.buildMenus();

            controls = uipanel('Parent', obj.Figure, 'Title', pd_ui_text('Data and run'), ...
                'Units', 'normalized', 'Position', [0.01, 0.02, 0.265, 0.96]);
            tabs = uitabgroup('Parent', obj.Figure, 'Units', 'normalized', ...
                'Position', [0.285, 0.02, 0.705, 0.96]);
            obj.TabGroup = tabs;
            plotTab = uitab('Parent', tabs, 'Title', pd_ui_text('Result plots'));
            obj.VisualizationTab = plotTab;
            plotDataTab = uitab('Parent', tabs, 'Title', pd_ui_text('Plot data'));
            obj.PlotDataTab = plotDataTab;
            computeTab = uitab('Parent', tabs, 'Title', pd_ui_text('Calculation parameters'));
            obj.ComputeTab = computeTab;
            plotOptionTab = uitab('Parent', tabs, 'Title', pd_ui_text('Plot settings'));
            outputTab = uitab('Parent', tabs, 'Title', pd_ui_text('Output settings'));
            resultTab = uitab('Parent', tabs, 'Title', pd_ui_text('Result summary'));

            addLabel(controls, pd_ui_text('Data file'), 0.925);
            obj.FileEdit = addEdit(controls, '', [0.04, 0.875, 0.55, 0.048]);
            set(obj.FileEdit, 'TooltipString', pd_ui_text('Select a chunk, cluster, or distribution data file'));
            addButton(controls, pd_ui_text('Browse'), [0.61, 0.875, 0.16, 0.048], ...
                @(src, evt)obj.browseFile(src, evt));
            sampleButton = addButton(controls, pd_ui_text('Example data'), [0.79, 0.875, 0.17, 0.048], ...
                @(src, evt)obj.loadDefaultExample(src, evt));
            set(sampleButton, 'TooltipString', pd_ui_text('Load ready-to-run test data for the current analysis type'));

            addLabel(controls, pd_ui_text('Analysis type'), 0.815);
            obj.TaskValues = {'chunk','network2d','cluster','vx','massx'};
            taskLabels = cellfun(@analysisDescription, obj.TaskValues, ...
                'UniformOutput', false);
            obj.TaskPopup = addPopup(controls, taskLabels, ...
                [0.04, 0.765, 0.92, 0.048], @(src, evt)obj.taskChanged(src, evt));

            addLabel(controls, pd_ui_text('Timestep selection mode'), 0.705);
            obj.SelectionPopup = addPopup(controls, {'Index','TimeStep','Time'}, ...
                [0.04, 0.655, 0.92, 0.048], @(src, evt)obj.selectionModeChanged(src, evt));

            addLabel(controls, pd_ui_text('Timestep/index'), 0.595);
            obj.SelectionEdit = addEdit(controls, '1', [0.04, 0.545, 0.92, 0.048]);

            addLabel(controls, pd_ui_text('Slurm file (required for physical-time selection)'), 0.485);
            obj.SlurmEdit = addEdit(controls, '', [0.04, 0.435, 0.72, 0.048]);
            obj.SlurmBrowseButton = addButton(controls, 'Browse', [0.78, 0.435, 0.18, 0.048], ...
                @(src, evt)obj.browseSlurm(src, evt));

            uicontrol('Parent', controls, 'Style', 'text', 'String', pd_ui_text('Progress display'), ...
                'Units', 'normalized', 'Position', [0.04, 0.375, 0.44, 0.035], ...
                'HorizontalAlignment', 'left');
            uicontrol('Parent', controls, 'Style', 'text', 'String', pd_ui_text('Result retention level'), ...
                'Units', 'normalized', 'Position', [0.52, 0.375, 0.44, 0.035], ...
                'HorizontalAlignment', 'left');
            obj.ProgressPopup = addPopup(controls, {'auto','console','waitbar','off'}, ...
                [0.04, 0.325, 0.44, 0.048], []);
            obj.ResultLevelPopup = addPopup(controls, {'standard','summary','full'}, ...
                [0.52, 0.325, 0.44, 0.048], []);

            obj.RunButton = addButton(controls, pd_ui_text('Run analysis (F5)'), [0.04, 0.255, 0.68, 0.058], ...
                @(src, evt)obj.runCurrentAnalysis(src, evt));
            set(obj.RunButton, 'FontWeight', 'bold');
            set(obj.RunButton, 'TooltipString', pd_ui_text('Run the current analysis and display result plots when complete'));
            obj.CancelButton = addButton(controls, pd_ui_text('Cancel (Esc)'), [0.74, 0.255, 0.22, 0.058], ...
                @(src, evt)obj.requestCancellation(src, evt));
            set(obj.CancelButton, 'Enable', 'off');
            addButton(controls, pd_ui_text('Save configuration'), [0.04, 0.195, 0.44, 0.045], ...
                @(src, evt)obj.saveConfiguration(src, evt));
            addButton(controls, pd_ui_text('Load configuration'), [0.52, 0.195, 0.44, 0.045], ...
                @(src, evt)obj.loadConfiguration(src, evt));
            addButton(controls, pd_ui_text('Restore default parameters'), [0.04, 0.14, 0.92, 0.045], ...
                @(src, evt)obj.resetAllOptions(src, evt));

            obj.StatusText = uicontrol('Parent', controls, 'Style', 'edit', ...
                'String', pd_ui_text('Ready'), 'Max', 3, 'Min', 0, 'Enable', 'inactive', ...
                'HorizontalAlignment', 'left', 'Units', 'normalized', ...
                'Position', [0.04, 0.025, 0.92, 0.10], 'BackgroundColor', 'white');

            uicontrol('Parent', plotTab, 'Style', 'text', 'String', pd_ui_text('Current view'), ...
                'Units', 'normalized', 'Position', [0.02, 0.935, 0.08, 0.035], ...
                'HorizontalAlignment', 'left');
            obj.PlotViewPopup = addPopup(plotTab, {pd_ui_text('All views')}, ...
                [0.10, 0.925, 0.19, 0.05], @(src, evt)obj.plotViewChanged(src, evt));
            obj.PlotViewIds = {'all'};
            obj.PlotViewCountText = uicontrol('Parent', plotTab, 'Style', 'text', ...
                'String', pd_ui_text('Not run yet'), 'Units', 'normalized', ...
                'Position', [0.30, 0.935, 0.16, 0.035]);
            obj.OpenSeparateFiguresCheckbox = uicontrol('Parent', plotTab, ...
                'Style', 'checkbox', 'String', pd_ui_text('Open separate figures after run'), ...
                'Value', 1, 'Units', 'normalized', ...
                'Position', [0.47, 0.925, 0.17, 0.05], ...
                'TooltipString', pd_ui_text('Open one scalable, savable MATLAB figure for each result by default'));
            addButton(plotTab, pd_ui_text('Open all figures'), [0.65, 0.925, 0.16, 0.05], ...
                @(src, evt)obj.openAllPlotFigures(src, evt));
            addButton(plotTab, pd_ui_text('Close separate figures'), [0.82, 0.925, 0.16, 0.05], ...
                @(src, evt)obj.closePlotFiguresFromButton(src, evt));
            obj.PlotFigureHandles = gobjects(1, 0);
            obj.PlotFigureAxes = gobjects(1, 0);
            obj.PlotFigureViewIds = {};
            obj.PlotPanel = uipanel('Parent', plotTab, 'BorderType', 'none', ...
                'Units', 'normalized', 'Position', [0.01, 0.11, 0.98, 0.80]);
            obj.PlotAxes = axes('Parent', obj.PlotPanel, 'Units', 'normalized', ...
                'Position', [0.09, 0.10, 0.84, 0.84]);
            addButton(plotTab, pd_ui_text('View plotted data'), [0.14, 0.035, 0.28, 0.055], ...
                @(src, evt)obj.showPlotDataTab(src, evt));
            addButton(plotTab, pd_ui_text('Apply plot settings and refresh'), [0.58, 0.035, 0.28, 0.055], ...
                @(src, evt)obj.replotCurrentResult(src, evt));

            obj.PlotDataPageSize = 500;
            obj.PlotDataPage = 1;
            obj.PlotDataSelectedRows = [];
            obj.CurrentPlotData = [];
            uicontrol('Parent', plotDataTab, 'Style', 'text', ...
                'String', pd_ui_text('Data view'), 'Units', 'normalized', ...
                'Position', [0.02, 0.935, 0.08, 0.035], ...
                'HorizontalAlignment', 'left');
            obj.PlotDataViewPopup = addPopup(plotDataTab, {pd_ui_text('Not run yet')}, ...
                [0.10, 0.925, 0.27, 0.05], ...
                @(src, evt)obj.plotDataViewChanged(src, evt));
            obj.PlotDataViewIds = {};
            obj.PlotDataPreviousButton = addButton(plotDataTab, ...
                pd_ui_text('Previous page'), [0.39, 0.925, 0.12, 0.05], ...
                @(src, evt)obj.changePlotDataPage(-1));
            obj.PlotDataPageText = uicontrol('Parent', plotDataTab, ...
                'Style', 'text', 'String', pd_ui_text('Page 0 of 0'), ...
                'Units', 'normalized', 'Position', [0.52, 0.935, 0.12, 0.035]);
            obj.PlotDataNextButton = addButton(plotDataTab, ...
                pd_ui_text('Next page'), [0.65, 0.925, 0.12, 0.05], ...
                @(src, evt)obj.changePlotDataPage(1));
            obj.PlotDataCopyHeadersCheckbox = uicontrol('Parent', plotDataTab, ...
                'Style', 'checkbox', ...
                'String', pd_ui_text('Include column headers when copying'), ...
                'Value', 1, 'Units', 'normalized', ...
                'Position', [0.79, 0.925, 0.19, 0.05], ...
                'TooltipString', pd_ui_text( ...
                    'Clear this option to copy numeric rows without column names.'));
            obj.PlotDataInfoText = uicontrol('Parent', plotDataTab, ...
                'Style', 'text', 'String', pd_ui_text('Run an analysis to inspect plot data.'), ...
                'Units', 'normalized', 'Position', [0.02, 0.875, 0.96, 0.04], ...
                'HorizontalAlignment', 'left');
            obj.PlotDataTable = uitable('Parent', plotDataTab, ...
                'Units', 'normalized', 'Position', [0.02, 0.12, 0.96, 0.74], ...
                'ColumnName', {}, 'ColumnEditable', false, ...
                'Data', cell(0, 0), ...
                'CellSelectionCallback', @(src, evt)obj.selectPlotDataRows(src, evt));
            addButton(plotDataTab, pd_ui_text('Copy selected rows'), ...
                [0.10, 0.035, 0.22, 0.055], ...
                @(src, evt)obj.copyPlotData(true));
            addButton(plotDataTab, pd_ui_text('Copy all data'), ...
                [0.39, 0.035, 0.22, 0.055], ...
                @(src, evt)obj.copyPlotData(false));
            addButton(plotDataTab, pd_ui_text('Export this view as CSV'), ...
                [0.68, 0.035, 0.22, 0.055], ...
                @(src, evt)obj.exportCurrentPlotData(src, evt));
            set([obj.PlotDataPreviousButton,obj.PlotDataNextButton], ...
                'Enable', 'off');

            obj.ComputeHeader = uicontrol('Parent', computeTab, 'Style', 'text', ...
                'String', '', 'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
                'Units', 'normalized', 'Position', [0.02, 0.94, 0.96, 0.04]);
            obj.ComputeTable = makeOptionTable(computeTab);
            set(obj.ComputeTable, 'Position', [0.02, 0.31, 0.96, 0.62], ...
                'CellSelectionCallback', @(src, evt)obj.selectComputeParameter(src, evt));
            editorPanel = uipanel('Parent', computeTab, 'Title', pd_ui_text('Parameter editor (select a table row to edit it here)'), ...
                'Units', 'normalized', 'Position', [0.02, 0.02, 0.96, 0.27]);
            obj.ParameterNameText = uicontrol('Parent', editorPanel, 'Style', 'text', ...
                'String', pd_ui_text('Select a calculation parameter first'), 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', 'Units', 'normalized', ...
                'Position', [0.02, 0.72, 0.42, 0.20]);
            obj.ParameterHelpText = uicontrol('Parent', editorPanel, 'Style', 'text', ...
                'String', '', 'HorizontalAlignment', 'left', 'Units', 'normalized', ...
                'Position', [0.02, 0.40, 0.94, 0.28]);
            obj.ParameterEdit = addEdit(editorPanel, '', [0.02, 0.10, 0.62, 0.24]);
            obj.ParameterPopup = addPopup(editorPanel, {' '}, [0.02, 0.10, 0.62, 0.24], []);
            obj.ParameterCheckbox = uicontrol('Parent', editorPanel, 'Style', 'checkbox', ...
                'String', pd_ui_text('Enabled'), 'Units', 'normalized', 'Position', [0.02, 0.10, 0.62, 0.24]);
            obj.ParameterApplyButton = addButton(editorPanel, pd_ui_text('Apply parameter value'), ...
                [0.70, 0.10, 0.27, 0.24], @(src, evt)obj.applySelectedParameter(src, evt));
            set([obj.ParameterEdit,obj.ParameterPopup,obj.ParameterCheckbox], 'Visible', 'off');
            set(obj.ParameterApplyButton, 'Enable', 'off');
            obj.PlotTable = makeOptionTable(plotOptionTab);
            uicontrol('Parent', plotOptionTab, 'Style', 'text', ...
                'String', pd_ui_text('Publication style preset'), 'HorizontalAlignment', 'left', ...
                'Units', 'normalized', 'Position', [0.02, 0.945, 0.13, 0.035]);
            obj.PlotPresetPopup = addPopup(plotOptionTab, ...
                {pd_ui_text('Journal color'),pd_ui_text('Monochrome print'),pd_ui_text('High contrast'),pd_ui_text('Restore defaults')}, ...
                [0.15, 0.935, 0.20, 0.05], []);
            addButton(plotOptionTab, pd_ui_text('Apply preset'), [0.37, 0.935, 0.14, 0.05], ...
                @(src, evt)obj.applyPlotPreset(src, evt));
            uicontrol('Parent', plotOptionTab, 'Style', 'text', ...
                'String', pd_ui_text('All parameters can be refined in the table or editor below'), ...
                'HorizontalAlignment', 'left', 'Units', 'normalized', ...
                'Position', [0.54, 0.945, 0.43, 0.035]);
            set(obj.PlotTable, 'Position', [0.02, 0.31, 0.96, 0.60], ...
                'CellSelectionCallback', @(src, evt)obj.selectPlotParameter(src, evt));
            plotEditorPanel = uipanel('Parent', plotOptionTab, ...
                'Title', pd_ui_text('Plot settings editor (changes do not require recalculation)'), 'Units', 'normalized', ...
                'Position', [0.02, 0.02, 0.96, 0.27]);
            obj.PlotParameterNameText = uicontrol('Parent', plotEditorPanel, ...
                'Style', 'text', 'String', pd_ui_text('Select a plot setting first'), ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
                'Units', 'normalized', 'Position', [0.02, 0.72, 0.42, 0.20]);
            obj.PlotParameterHelpText = uicontrol('Parent', plotEditorPanel, ...
                'Style', 'text', 'String', '', 'HorizontalAlignment', 'left', ...
                'Units', 'normalized', 'Position', [0.02, 0.40, 0.94, 0.28]);
            obj.PlotParameterEdit = addEdit(plotEditorPanel, '', ...
                [0.02, 0.10, 0.62, 0.24]);
            obj.PlotParameterPopup = addPopup(plotEditorPanel, {' '}, ...
                [0.02, 0.10, 0.62, 0.24], []);
            obj.PlotParameterCheckbox = uicontrol('Parent', plotEditorPanel, ...
                'Style', 'checkbox', 'String', pd_ui_text('Enabled'), 'Units', 'normalized', ...
                'Position', [0.02, 0.10, 0.62, 0.24]);
            obj.PlotParameterApplyButton = addButton(plotEditorPanel, ...
                pd_ui_text('Apply plot value'), [0.70, 0.10, 0.27, 0.24], ...
                @(src, evt)obj.applySelectedPlotParameter(src, evt));
            set([obj.PlotParameterEdit,obj.PlotParameterPopup, ...
                obj.PlotParameterCheckbox], 'Visible', 'off');
            set(obj.PlotParameterApplyButton, 'Enable', 'off');
            obj.OutputTable = makeOptionTable(outputTab);
            set(obj.ComputeTable, 'CellEditCallback', @(src, evt)obj.validateEditableTables(src, evt));
            set(obj.PlotTable, 'CellEditCallback', @(src, evt)obj.plotTableEdited(src, evt));
            set(obj.OutputTable, 'CellEditCallback', @(src, evt)obj.validateEditableTables(src, evt));
            addButton(outputTab, pd_ui_text('Select output directory'), [0.16, 0.035, 0.25, 0.055], ...
                @(src, evt)obj.browseOutputDirectory(src, evt));
            addButton(outputTab, pd_ui_text('Export current result'), [0.58, 0.035, 0.25, 0.055], ...
                @(src, evt)obj.exportCurrentResult(src, evt));

            obj.ResultTable = uitable('Parent', resultTab, 'Units', 'normalized', ...
                'Position', [0.02, 0.03, 0.96, 0.94], ...
                'ColumnName', {pd_ui_text('Result field'),pd_ui_text('Value')}, ...
                'ColumnEditable', [false false], 'ColumnWidth', {360 260}, ...
                'Data', cell(0, 2));
        end

        function buildMenus(obj)
            sampleMenu = uimenu(obj.Figure, 'Label', pd_ui_text('Example data'));
            uimenu(sampleMenu, 'Label', pd_ui_text('1D field (401 grid points x 3 timesteps)'), ...
                'Callback', @(src, evt)obj.loadExample('chunk1d'));
            uimenu(sampleMenu, 'Label', pd_ui_text('2D field (61 x 41 grid x 3 timesteps)'), ...
                'Callback', @(src, evt)obj.loadExample('chunk2d'));
            uimenu(sampleMenu, 'Label', pd_ui_text('Particle/cluster (1200 items x 3 timesteps)'), ...
                'Callback', @(src, evt)obj.loadExample('cluster'));
            uimenu(sampleMenu, 'Label', pd_ui_text('mass-x cumulative distribution'), ...
                'Callback', @(src, evt)obj.loadExample('massx'), 'Separator', 'on');
            uimenu(sampleMenu, 'Label', pd_ui_text('mass-x 2D sliced distribution'), ...
                'Callback', @(src, evt)obj.loadExample('massx2d'));
            uimenu(sampleMenu, 'Label', pd_ui_text('2D pore network'), ...
                'Callback', @(src, evt)obj.loadExample('network2d'));
            uimenu(sampleMenu, 'Label', pd_ui_text('mass-v cumulative distribution'), ...
                'Callback', @(src, evt)obj.loadExample('vx'));

            toolsMenu = uimenu(obj.Figure, 'Label', pd_ui_text('Tools'));
            uimenu(toolsMenu, 'Label', pd_ui_text('Validate current input and parameters'), ...
                'Callback', @(src, evt)obj.preflightCurrentInput(src, evt));
            uimenu(toolsMenu, 'Label', pd_ui_text('System diagnostics'), ...
                'Callback', @(src, evt)obj.showSystemDiagnostics(src, evt));
            uimenu(toolsMenu, 'Label', pd_ui_text('Open run log'), 'Separator', 'on', ...
                'Callback', @(src, evt)obj.openRunLog(src, evt));

            helpMenu = uimenu(obj.Figure, 'Label', pd_ui_text('Help'));
            uimenu(helpMenu, 'Label', pd_ui_text('Quick start and shortcuts'), ...
                'Callback', @(src, evt)obj.showQuickStart(src, evt));
            uimenu(helpMenu, 'Label', pd_ui_text('About POST_DATA2'), ...
                'Callback', @(src, evt)obj.showAbout(src, evt));
        end

        function loadDefaultExample(obj, ~, ~)
            task = obj.currentAnalysisType();
            if strcmp(task, 'chunk')
                data = get(obj.ComputeTable, 'Data');
                row = find(strcmp('ChunkDim', data(:, 1)), 1, 'first');
                if ~isempty(row) && strcmpi(strtrim(pd_to_char(data{row, 2})), '1d')
                    kind = 'chunk1d';
                else
                    kind = 'chunk2d';
                end
            else
                kind = task;
            end
            obj.loadExample(kind);
        end

        function task = currentAnalysisType(obj)
            index = get(obj.TaskPopup, 'Value');
            index = max(1, min(index, numel(obj.TaskValues)));
            task = obj.TaskValues{index};
        end

        function setAnalysisType(obj, value)
            value = pd_normalize_analysis_type(value);
            index = find(strcmp(value, obj.TaskValues), 1, 'first');
            if isempty(index)
                error('postdata:BadPopupValue', ...
                    'Unsupported analysis type: %s', value);
            end
            set(obj.TaskPopup, 'Value', index);
        end

        function loadExample(obj, kind)
            rootDir = fileparts(which('postdata_app'));
            generatedDir = fullfile(rootDir, 'fixtures', 'generated');
            sampleDir = fullfile(rootDir, 'fixtures', 'sample');
            switch lower(kind)
                case 'chunk1d'
                    task = 'chunk';
                    filePath = fullfile(generatedDir, 'large_bin1d_dx_0.025.txt');
                    overrides = {'ChunkDim', '1d'};
                case 'chunk2d'
                    task = 'chunk';
                    filePath = fullfile(generatedDir, ...
                        'large_bin2d_dx_0.05_dy_0.06_Lz_1.txt');
                    overrides = {'ChunkDim', '2d'};
                case 'cluster'
                    task = 'cluster';
                    filePath = fullfile(generatedDir, 'large_cluster_chunk.txt');
                    overrides = {'Dim', 2, 'Dx', 0.025, 'MeanNumBins', 30};
                case 'massx'
                    task = 'massx';
                    filePath = fullfile(generatedDir, 'large_bin1d_dx_0.025.txt');
                    overrides = {'ChunkDim', '1d', 'SphDimension', 2, ...
                        'InitialDensity', 7.3, ...
                        'ParticleSpacing', 0.005, 'TransverseWidth', 2.46, ...
                        'RawLengthUnitUm', 10, 'CoordinateFactor', 10};
                case 'massx2d'
                    task = 'massx';
                    filePath = fullfile(generatedDir, ...
                        'large_bin2d_dx_0.05_dy_0.06_Lz_1.txt');
                    overrides = {'ChunkDim', '2d', 'SphDimension', 2, ...
                        'InitialDensity', 7.3, ...
                        'ParticleSpacing', 0.01, 'RawLengthUnitUm', 10, ...
                        'CoordinateFactor', 10, ...
                        'SliceCentersY', [0.6 1.2 1.8], ...
                        'SliceWidthsY', 0.6};
                case 'network2d'
                    task = 'network2d';
                    filePath = fullfile(generatedDir, ...
                        'large_bin2d_dx_0.05_dy_0.06_Lz_1.txt');
                    overrides = {'GeometryMode', 'original', ...
                        'PositionAxis', 'both', 'ProfileAxis', 'both'};
                case 'vx'
                    task = 'vx';
                    filePath = fullfile(sampleDir, 'vx_chunk_test.txt');
                    overrides = {};
                otherwise
                    error('postdata:UnknownExample', 'Unknown example: %s', kind);
            end
            if ~exist(filePath, 'file')
                error('postdata:MissingExampleData', ...
                    'Example data is missing: %s', filePath);
            end
            obj.setAnalysisType(task);
            obj.taskChanged([], []);
            for i = 1:2:numel(overrides)
                obj.setComputeOptionValue(overrides{i}, overrides{i + 1});
            end
            set(obj.FileEdit, 'String', filePath, 'BackgroundColor', 'white');
            setPopupValue(obj.SelectionPopup, 'Index');
            set(obj.SelectionEdit, 'String', '1');
            obj.selectionModeChanged([], []);
            set(obj.TabGroup, 'SelectedTab', obj.ComputeTab);
            obj.updateRunButtonForParameters();
            obj.setStatus(sprintf(pd_ui_text('Loaded %s example. Review parameters and press F5 to run.'), ...
                analysisDescription(task)), 'success');
        end

        function setComputeOptionValue(obj, name, value)
            data = get(obj.ComputeTable, 'Data');
            row = find(strcmp(name, data(:, 1)), 1, 'first');
            if isempty(row)
                error('postdata:MissingExampleOption', ...
                    'Example option is unavailable: %s', name);
            end
            data{row, 2} = pd_format_option_editor_value( ...
                value, obj.OptionCatalog(row).Type);
            set(obj.ComputeTable, 'Data', data);
            obj.updateDependencyDisplay();
        end

        function showQuickStart(~, ~, ~)
            helpdlg({pd_ui_text('Quick workflow:'), ...
                pd_ui_text('1. Select a data file, or load one from the Example Data menu.'), ...
                pd_ui_text('2. Select an analysis type and review the Calculation Parameters tab.'), ...
                pd_ui_text('3. Press F5 to run; the Result Plots tab opens automatically.'), ...
                pd_ui_text('4. A separate MATLAB figure opens for each result by default.'), ...
                pd_ui_text('5. Adjust and export from Plot Settings and Output Settings.'), ...
                pd_ui_text('Shortcuts: F5 run, Esc cancel, Ctrl+O browse, Ctrl+E load example.')}, ...
                pd_ui_text('POST_DATA2 Quick Start'));
        end

        function showAbout(~, ~, ~)
            info = pd_version();
            helpdlg({['POST_DATA2 ', info.version], ...
                pd_ui_text('MATLAB post-processing software'), ...
                [pd_ui_text('Minimum compatible release:'), info.minimumMatlabRelease], ...
                pd_ui_text('Supports chunk, cluster, mass-v, mass-x, and network2d.')}, ...
                pd_ui_text('About POST_DATA2'));
        end

        function preflightCurrentInput(obj, ~, ~)
            try
                request = pd_validate_request(obj.buildRequest());
                filePath = pd_resolve_input_file(request);
                report = pd_preflight_input(request, filePath);
                shown = report.variables(1:min(12, numel(report.variables)));
                variableText = strjoin(shown, ', ');
                if numel(report.variables) > numel(shown)
                    variableText = sprintf(pd_ui_text('%s, ... (%d columns total)'), ...
                        variableText, numel(report.variables));
                end
                helpdlg({ ...
                    pd_ui_text('Input and parameter validation passed'), ...
                    [pd_ui_text('File:'), report.filePath], ...
                    sprintf(pd_ui_text('Size: %.3f MB'), report.fileSize / 1024 / 1024), ...
                    sprintf(pd_ui_text('First timestep: %g'), report.firstTimestep), ...
                    sprintf(pd_ui_text('First-block data rows: %d'), report.firstBlockRows), ...
                    [pd_ui_text('Variables:'), variableText]}, pd_ui_text('POST_DATA2 Preflight Check'));
                obj.setStatus(pd_ui_text('Input file, analysis type, and parameters passed validation.'), 'success');
            catch err
                obj.showError(err);
            end
        end

        function showSystemDiagnostics(obj, ~, ~)
            try
                report = pd_system_diagnostics();
                if report.imageProcessingToolboxAvailable
                    imageToolboxText = pd_ui_text('Available');
                else
                    imageToolboxText = pd_ui_text('Unavailable (only optional skeleton/thickness metrics are affected)');
                end
                lines = { ...
                    sprintf('POST_DATA2 %s', report.software.version), ...
                    [pd_ui_text('MATLAB:'), report.matlabVersion], ...
                    [pd_ui_text('Runtime release:'), report.matlabRelease], ...
                    [pd_ui_text('Project directory:'), report.projectRoot], ...
                    [pd_ui_text('Default output:'), report.defaultOutputDirectory], ...
                    [pd_ui_text('Run log:'), report.logFile], ...
                    [pd_ui_text('Image Processing Toolbox:'), imageToolboxText]};
                if isempty(report.warnings)
                    lines{end + 1} = pd_ui_text('Conclusion: core runtime checks passed.');
                    state = 'success';
                else
                    lines{end + 1} = pd_ui_text('Notes:');
                    for i = 1:numel(report.warnings)
                        lines{end + 1} = ['- ', report.warnings{i}]; %#ok<AGROW>
                    end
                    state = 'normal';
                end
                helpdlg(lines, pd_ui_text('POST_DATA2 System Diagnostics'));
                obj.setStatus(pd_ui_text('System diagnostics completed.'), state);
            catch err
                obj.showError(err);
            end
        end

        function openRunLog(obj, ~, ~)
            filePath = pd_default_log_file();
            if ~exist(filePath, 'file')
                fid = fopen(filePath, 'a');
                if fid >= 0, fclose(fid); end
            end
            try
                open(filePath);
                obj.setStatus([pd_ui_text('Opened run log:'), filePath]);
            catch err
                obj.showError(err);
            end
        end

        function handleShortcut(obj, ~, event)
            modifiers = event.Modifier;
            if ischar(modifiers), modifiers = {modifiers}; end
            controlDown = any(strcmpi(modifiers, 'control'));
            if strcmpi(event.Key, 'f5') || ...
                    (controlDown && strcmpi(event.Key, 'r'))
                obj.runCurrentAnalysis([], []);
            elseif strcmpi(event.Key, 'escape')
                if strcmpi(char(get(obj.CancelButton, 'Enable')), 'on')
                    obj.requestCancellation([], []);
                end
            elseif controlDown && strcmpi(event.Key, 'o')
                obj.browseFile([], []);
            elseif controlDown && strcmpi(event.Key, 'e')
                obj.loadDefaultExample([], []);
            end
        end

        function browseFile(obj, ~, ~)
            initial = browseStart(get(obj.FileEdit, 'String'), '*.txt');
            [name, folder] = uigetfile( ...
                {'*.txt', 'POST_DATA text files (*.txt)'; '*.*', 'All files'}, ...
                pd_ui_text('Select POST_DATA data file'), initial);
            if isequal(name, 0), return; end
            set(obj.FileEdit, 'String', fullfile(folder, name), 'BackgroundColor', 'white');
            obj.setStatus(pd_ui_text('Data file selected; review the analysis type and parameters.'));
        end

        function browseSlurm(obj, ~, ~)
            initial = browseStart(get(obj.SlurmEdit, 'String'), '*.log');
            [name, folder] = uigetfile( ...
                {'*.log;*.txt', 'Slurm files'; '*.*', 'All files'}, ...
                pd_ui_text('Select Slurm log'), initial);
            if isequal(name, 0), return; end
            set(obj.SlurmEdit, 'String', fullfile(folder, name));
        end

        function browseOutputDirectory(obj, ~, ~)
            current = obj.getOutputOptions();
            folder = uigetdir(current.Directory, 'Choose output directory');
            if isequal(folder, 0), return; end
            data = get(obj.OutputTable, 'Data');
            row = find(strcmp('Directory', data(:, 1)), 1, 'first');
            data{row, 2} = folder;
            set(obj.OutputTable, 'Data', data);
        end

        function taskChanged(obj, source, ~)
            task = obj.currentAnalysisType();
            if ~isempty(obj.ActiveAnalysisType) && ishghandle(obj.ComputeTable)
                obj.OptionDataCache.(obj.ActiveAnalysisType) = get(obj.ComputeTable, 'Data');
            end
            obj.OptionCatalog = pd_option_catalog(task);
            if isfield(obj.OptionDataCache, task) && ...
                    size(obj.OptionDataCache.(task), 1) == numel(obj.OptionCatalog)
                data = obj.OptionDataCache.(task);
            else
                data = pd_catalog_table_data(obj.OptionCatalog);
            end
            data = pd_normalize_catalog_table_data(obj.OptionCatalog, data);
            set(obj.ComputeTable, 'Data', data);
            obj.ActiveAnalysisType = task;
            obj.SelectedParameterIndex = [];
            obj.updateDependencyDisplay();
            set(obj.ComputeHeader, 'String', computeHeaderText( ...
                task, numel(obj.OptionCatalog)));
            if strcmp(task, 'massx')
                row = find(strcmp('InitialDensity', {obj.OptionCatalog.Name}), ...
                    1, 'first');
                if ~isempty(row)
                    obj.SelectedParameterIndex = row;
                    obj.showSelectedParameter();
                end
            end
            if ~isempty(source) && ishghandle(obj.TabGroup) && ishghandle(obj.ComputeTab)
                set(obj.TabGroup, 'SelectedTab', obj.ComputeTab);
            end
            obj.updateRunButtonForParameters();
            obj.setStatus(sprintf(pd_ui_text('Switched to %s with %d editable calculation parameters.'), ...
                analysisDescription(task), numel(obj.OptionCatalog)));
        end

        function resetAllOptions(obj, ~, ~)
            obj.OptionDataCache = struct();
            obj.ActiveAnalysisType = '';
            obj.taskChanged([], []);
            obj.PlotCatalog = pd_plot_option_catalog();
            set(obj.PlotTable, 'Data', pd_catalog_table_data(obj.PlotCatalog));
            obj.SelectedPlotParameterIndex = [];
            obj.showSelectedPlotParameter();
            obj.OutputCatalog = pd_output_option_catalog();
            set(obj.OutputTable, 'Data', pd_catalog_table_data(obj.OutputCatalog));
            obj.selectionModeChanged([], []);
        end

        function runCurrentAnalysis(obj, ~, ~)
            obj.CancelRequested = false;
            set(obj.RunButton, 'Enable', 'off');
            set(obj.CancelButton, 'Enable', 'on');
            cleanupObj = onCleanup(@() obj.finishExecutionState());
            try
                request = obj.buildRequest();
                obj.setStatus(pd_ui_text('Reading data and running analysis, please wait...'), 'busy');
                drawnow;
                obj.CurrentResult = postdata_run(request);
                obj.updatePlotViewChoices();
                obj.replotCurrentResult([], []);
                set(obj.TabGroup, 'SelectedTab', obj.VisualizationTab);
                summary = pd_result_summary(obj.CurrentResult);
                set(obj.ResultTable, 'Data', summary);
                obj.setStatus(sprintf(pd_ui_text('%s completed: timestep = %g, generated %d result views.'), ...
                    analysisDescription(obj.CurrentResult.analysisType), ...
                    obj.CurrentResult.timestep, numel(pd_result_plot_views(obj.CurrentResult))), ...
                    'success');
                output = obj.getOutputOptions();
                if output.AutoExport
                    obj.exportWithOptions(output);
                end
            catch err
                if strcmp(err.identifier, 'postdata:UserCancelled')
                    obj.setStatus(pd_ui_text('Analysis cancelled by user.'));
                else
                    obj.showError(err);
                end
            end
        end

        function request = buildRequest(obj)
            filePath = strtrim(get(obj.FileEdit, 'String'));
            if isempty(filePath) || ~exist(filePath, 'file')
                set(obj.FileEdit, 'BackgroundColor', [1.0, 0.84, 0.84]);
                error('postdata:MissingInputFile', pd_ui_text('Select an existing data file.'));
            end
            set(obj.FileEdit, 'BackgroundColor', 'white');
            selectionValue = str2double(get(obj.SelectionEdit, 'String'));
            if ~isfinite(selectionValue)
                error('postdata:BadSelectionValue', 'Selection value must be numeric.');
            end
            task = obj.currentAnalysisType();
            if ~strcmp(task, obj.ActiveAnalysisType)
                obj.taskChanged([], []);
            end
            request = pd_create_request(task);
            request.filePath = filePath;
            request.baseDir = fileparts(filePath);
            request.selection.mode = popupValue(obj.SelectionPopup);
            request.selection.value = selectionValue;
            request.selection.slurmPath = strtrim(get(obj.SlurmEdit, 'String'));
            request.progressMode = popupValue(obj.ProgressPopup);
            request.resultLevel = popupValue(obj.ResultLevelPopup);
            request.makePlots = false;
            request.plotOptions = pd_plot_options_from_table(obj.PlotCatalog, ...
                get(obj.PlotTable, 'Data'));
            request.execution.cancelCallback = @() obj.isCancellationRequested();
            request.execution.progressCallback = @(fraction, message) ...
                obj.updateExecutionProgress(fraction, message);
            request = pd_apply_analysis_options(request, obj.OptionCatalog, ...
                get(obj.ComputeTable, 'Data'));
        end

        function requestCancellation(obj, ~, ~)
            obj.CancelRequested = true;
            set(obj.CancelButton, 'Enable', 'off');
            obj.setStatus('Cancellation requested; waiting for the current checkpoint...');
            drawnow;
        end

        function value = isCancellationRequested(obj)
            value = ~isempty(obj.CancelRequested) && obj.CancelRequested;
        end

        function updateExecutionProgress(obj, fraction, message)
            obj.setStatus(sprintf('%s (%.0f%%)', message, 100 * fraction));
            drawnow;
        end

        function finishExecutionState(obj)
            obj.updateRunButtonForParameters();
            if ishghandle(obj.CancelButton), set(obj.CancelButton, 'Enable', 'off'); end
            obj.CancelRequested = false;
        end

        function updateRunButtonForParameters(obj)
            state = 'off';
            try
                if ~isempty(obj.OptionCatalog) && ishghandle(obj.ComputeTable)
                    pd_validate_catalog_values(obj.OptionCatalog, ...
                        get(obj.ComputeTable, 'Data'));
                    state = 'on';
                end
            catch
                state = 'off';
            end
            if ishghandle(obj.RunButton)
                set(obj.RunButton, 'Enable', state);
            end
        end

        function validateEditableTables(obj, source, event)
            try
                if isequal(source, obj.ComputeTable) && ~isempty(event) && ~isempty(event.Indices)
                    row = event.Indices(1);
                    enabled = pd_catalog_enabled_mask(obj.OptionCatalog, get(obj.ComputeTable, 'Data'));
                    if ~enabled(row)
                        data = get(obj.ComputeTable, 'Data');
                        data{row, 2} = event.PreviousData;
                        set(obj.ComputeTable, 'Data', data);
                        obj.updateDependencyDisplay();
                        obj.setStatus('This parameter is inactive because its dependency is not enabled.');
                        return;
                    end
                end
                task = obj.currentAnalysisType();
                request = pd_create_request(task);
                pd_apply_analysis_options(request, obj.OptionCatalog, get(obj.ComputeTable, 'Data'));
                pd_plot_options_from_table(obj.PlotCatalog, get(obj.PlotTable, 'Data'));
                pd_output_options_from_table(obj.OutputCatalog, get(obj.OutputTable, 'Data'));
                computeData = pd_normalize_catalog_table_data( ...
                    obj.OptionCatalog, get(obj.ComputeTable, 'Data'));
                plotData = pd_normalize_catalog_table_data( ...
                    obj.PlotCatalog, get(obj.PlotTable, 'Data'));
                outputData = pd_normalize_catalog_table_data( ...
                    obj.OutputCatalog, get(obj.OutputTable, 'Data'));
                set(obj.ComputeTable, 'Data', computeData);
                set(obj.PlotTable, 'Data', plotData);
                set(obj.OutputTable, 'Data', outputData);
                if ishghandle(obj.RunButton), set(obj.RunButton, 'Enable', 'on'); end
                obj.setStatus('Parameter syntax is valid.');
                obj.updateDependencyDisplay();
            catch err
                if ishghandle(obj.RunButton), set(obj.RunButton, 'Enable', 'off'); end
                obj.setStatus([pd_ui_text('Invalid parameter:'), ' ', ...
                    localizeValidationMessage(err.message)]);
            end
        end

        function selectionModeChanged(obj, ~, ~)
            enabled = strcmpi(popupValue(obj.SelectionPopup), 'Time');
            if enabled, state = 'on'; else, state = 'off'; end
            set(obj.SlurmEdit, 'Enable', state);
            set(obj.SlurmBrowseButton, 'Enable', state);
        end

        function selectComputeParameter(obj, ~, event)
            if isempty(event.Indices), return; end
            obj.SelectedParameterIndex = event.Indices(1);
            obj.showSelectedParameter();
        end

        function showSelectedParameter(obj)
            set([obj.ParameterEdit,obj.ParameterPopup,obj.ParameterCheckbox], 'Visible', 'off');
            if isempty(obj.SelectedParameterIndex) || ...
                    obj.SelectedParameterIndex > numel(obj.OptionCatalog)
                set(obj.ParameterApplyButton, 'Enable', 'off');
                return;
            end
            row = obj.SelectedParameterIndex;
            meta = obj.OptionCatalog(row);
            data = get(obj.ComputeTable, 'Data');
            enabled = pd_catalog_enabled_mask(obj.OptionCatalog, data);
            set(obj.ParameterNameText, 'String', sprintf('%s (%s)', meta.Name, meta.Type));
            help = pd_ui_text(meta.Description);
            if ~enabled(row)
                help = sprintf(pd_ui_text('%s | Inactive until %s = %s'), help, ...
                    meta.DependsOn, pd_format_option_value(meta.DependsValue));
            end
            set(obj.ParameterHelpText, 'String', help);
            if strcmp(meta.Type, 'logical')
                value = pd_parse_option_value(data{row, 2}, meta.Type);
                set(obj.ParameterCheckbox, 'Value', double(value), 'Visible', 'on');
            elseif ~isempty(meta.AllowedValues)
                labels = meta.AllowedValues;
                labels(cellfun(@isempty, labels)) = {'<empty>'};
                set(obj.ParameterPopup, 'String', labels, 'Visible', 'on');
                current = data{row, 2};
                if ~ischar(current), current = pd_format_option_value(current); end
                index = find(strcmpi(current, meta.AllowedValues), 1, 'first');
                if isempty(index), index = 1; end
                set(obj.ParameterPopup, 'Value', index);
            else
                set(obj.ParameterEdit, 'String', data{row, 2}, 'Visible', 'on');
            end
            if enabled(row), state = 'on'; else, state = 'off'; end
            set(obj.ParameterApplyButton, 'Enable', state);
        end

        function applySelectedParameter(obj, ~, ~)
            if isempty(obj.SelectedParameterIndex), return; end
            row = obj.SelectedParameterIndex;
            meta = obj.OptionCatalog(row);
            if strcmp(meta.Type, 'logical')
                if get(obj.ParameterCheckbox, 'Value'), textValue = 'true'; else, textValue = 'false'; end
            elseif ~isempty(meta.AllowedValues)
                index = get(obj.ParameterPopup, 'Value');
                textValue = meta.AllowedValues{index};
            else
                textValue = get(obj.ParameterEdit, 'String');
            end
            data = get(obj.ComputeTable, 'Data');
            previous = data{row, 2};
            data{row, 2} = textValue;
            try
                pd_validate_catalog_values(obj.OptionCatalog, data);
                data = pd_normalize_catalog_table_data(obj.OptionCatalog, data);
                set(obj.ComputeTable, 'Data', data);
                obj.updateDependencyDisplay();
                obj.showSelectedParameter();
                obj.validateEditableTables([], []);
            catch err
                data{row, 2} = previous;
                set(obj.ComputeTable, 'Data', data);
                obj.showError(err);
            end
        end

        function selectPlotParameter(obj, ~, event)
            if isempty(event.Indices), return; end
            obj.SelectedPlotParameterIndex = event.Indices(1);
            obj.showSelectedPlotParameter();
        end

        function applyPlotPreset(obj, ~, ~)
            data = get(obj.PlotTable, 'Data');
            preset = popupValue(obj.PlotPresetPopup);
            switch preset
                case pd_ui_text('Journal color')
                    values = {'ColorPalette','colorblind'; ...
                        'SeriesStyleMode','color-and-style'; 'FontName','Times New Roman'; ...
                        'FontSize',11; 'TitleFontSize',12; 'TitleWeight','normal'; ...
                        'LineWidth',1.8; 'AxisLineWidth',1.0; 'MarkerSymbol','none'; ...
                        'ShowGrid',false; 'ShowMinorGrid',false; 'BoxOn',true; ...
                        'TickDirection','out'; 'LegendBox',false};
                case pd_ui_text('Monochrome print')
                    values = {'ColorPalette','grayscale'; ...
                        'SeriesStyleMode','monochrome'; 'FontName','Times New Roman'; ...
                        'FontSize',11; 'TitleFontSize',12; 'LineWidth',1.8; ...
                        'MarkerSymbol','none'; 'ShowGrid',false; 'BoxOn',true; ...
                        'TickDirection','out'; 'LegendBox',false};
                case pd_ui_text('High contrast')
                    values = {'ColorPalette','highcontrast'; ...
                        'SeriesStyleMode','color-and-style'; 'FontName','Arial'; ...
                        'FontSize',14; 'TitleFontSize',16; 'TitleWeight','bold'; ...
                        'LineWidth',2.4; 'MarkerSymbol','o'; 'MarkerSize',7; ...
                        'ShowGrid',true; 'GridAlpha',0.22};
                otherwise
                    data = pd_catalog_table_data(obj.PlotCatalog);
                    set(obj.PlotTable, 'Data', data);
                    obj.refreshPlotAppearance();
                    return;
            end
            for i = 1:size(values, 1)
                data = replaceTableOption(data, values{i, 1}, values{i, 2});
            end
            pd_plot_options_from_table(obj.PlotCatalog, data);
            set(obj.PlotTable, 'Data', data);
            obj.showSelectedPlotParameter();
            obj.refreshPlotAppearance();
        end

        function showSelectedPlotParameter(obj)
            set([obj.PlotParameterEdit,obj.PlotParameterPopup, ...
                obj.PlotParameterCheckbox], 'Visible', 'off');
            if isempty(obj.SelectedPlotParameterIndex) || ...
                    obj.SelectedPlotParameterIndex > numel(obj.PlotCatalog)
                set(obj.PlotParameterApplyButton, 'Enable', 'off');
                return;
            end
            row = obj.SelectedPlotParameterIndex;
            meta = obj.PlotCatalog(row);
            data = get(obj.PlotTable, 'Data');
            set(obj.PlotParameterNameText, 'String', ...
                sprintf('%s (%s)', meta.Name, meta.Type));
            set(obj.PlotParameterHelpText, 'String', meta.Description);
            if strcmp(meta.Type, 'logical')
                value = pd_parse_option_value(data{row, 2}, meta.Type);
                set(obj.PlotParameterCheckbox, 'Value', double(value), ...
                    'Visible', 'on');
            elseif ~isempty(meta.AllowedValues)
                set(obj.PlotParameterPopup, 'String', meta.AllowedValues, ...
                    'Visible', 'on');
                current = data{row, 2};
                if ~ischar(current), current = pd_format_option_value(current); end
                index = find(strcmpi(current, meta.AllowedValues), 1, 'first');
                if isempty(index), index = 1; end
                set(obj.PlotParameterPopup, 'Value', index);
            else
                set(obj.PlotParameterEdit, 'String', data{row, 2}, 'Visible', 'on');
            end
            set(obj.PlotParameterApplyButton, 'Enable', 'on');
        end

        function applySelectedPlotParameter(obj, ~, ~)
            if isempty(obj.SelectedPlotParameterIndex), return; end
            row = obj.SelectedPlotParameterIndex;
            meta = obj.PlotCatalog(row);
            if strcmp(meta.Type, 'logical')
                if get(obj.PlotParameterCheckbox, 'Value')
                    textValue = 'true';
                else
                    textValue = 'false';
                end
            elseif ~isempty(meta.AllowedValues)
                index = get(obj.PlotParameterPopup, 'Value');
                textValue = meta.AllowedValues{index};
            else
                textValue = get(obj.PlotParameterEdit, 'String');
            end
            data = get(obj.PlotTable, 'Data');
            previous = data{row, 2};
            data{row, 2} = textValue;
            set(obj.PlotTable, 'Data', data);
            try
                pd_plot_options_from_table(obj.PlotCatalog, data);
                obj.showSelectedPlotParameter();
                obj.validateEditableTables([], []);
                obj.refreshPlotAppearance();
            catch err
                data{row, 2} = previous;
                set(obj.PlotTable, 'Data', data);
                obj.showError(err);
            end
        end

        function plotTableEdited(obj, source, event)
            try
                pd_plot_options_from_table(obj.PlotCatalog, get(obj.PlotTable, 'Data'));
                obj.validateEditableTables(source, event);
                obj.refreshPlotAppearance();
            catch err
                if ~isempty(event) && ~isempty(event.Indices)
                    data = get(obj.PlotTable, 'Data');
                    data{event.Indices(1), event.Indices(2)} = event.PreviousData;
                    set(obj.PlotTable, 'Data', data);
                end
                obj.showError(err);
            end
        end

        function refreshPlotAppearance(obj)
            if isempty(obj.CurrentResult), return; end
            options = pd_plot_options_from_table(obj.PlotCatalog, ...
                get(obj.PlotTable, 'Data'));
            validAxes = obj.PlotAxes(ishghandle(obj.PlotAxes));
            if ~isempty(validAxes), pd_apply_publication_style(validAxes, options); end
            validSeparate = obj.PlotFigureAxes(ishghandle(obj.PlotFigureAxes));
            if ~isempty(validSeparate)
                pd_apply_publication_style(validSeparate, options);
            end
            drawnow;
            obj.setStatus(pd_ui_text('Plot settings applied immediately without recalculation.'), 'success');
        end

        function updateDependencyDisplay(obj)
            if isempty(obj.OptionCatalog), return; end
            data = get(obj.ComputeTable, 'Data');
            if size(data, 1) ~= numel(obj.OptionCatalog), return; end
            enabled = pd_catalog_enabled_mask(obj.OptionCatalog, data);
            for i = 1:numel(obj.OptionCatalog)
                description = pd_ui_text(obj.OptionCatalog(i).Description);
                if ~enabled(i)
                    description = sprintf(pd_ui_text('[inactive: %s=%s] %s'), ...
                        obj.OptionCatalog(i).DependsOn, ...
                        pd_format_option_value(obj.OptionCatalog(i).DependsValue), description);
                end
                data{i, 3} = description;
            end
            set(obj.ComputeTable, 'Data', data);
            if ~isempty(obj.SelectedParameterIndex), obj.showSelectedParameter(); end
        end

        function replotCurrentResult(obj, ~, ~)
            if isempty(obj.CurrentResult)
                return;
            end
            try
                options = pd_plot_options_from_table(obj.PlotCatalog, get(obj.PlotTable, 'Data'));
                popupIndex = get(obj.PlotViewPopup, 'Value');
                viewId = obj.PlotViewIds{popupIndex};
                views = pd_result_plot_views(obj.CurrentResult);
                if strcmp(viewId, 'all')
                    viewIds = {views.Id};
                else
                    viewIds = {viewId};
                end
                obj.preparePlotAxes(viewIds, options.UpdateMode);
                for i = 1:numel(viewIds)
                    pd_render_result(obj.PlotAxes(i), obj.CurrentResult, options, viewIds{i});
                end
                if get(obj.OpenSeparateFiguresCheckbox, 'Value') || ...
                        obj.hasAnySeparateFigures()
                    obj.renderSeparateFigures(options, false);
                end
                drawnow;
            catch err
                obj.showError(err);
            end
        end

        function plotViewChanged(obj, ~, ~)
            obj.replotCurrentResult([], []);
            if isempty(obj.CurrentResult) || isempty(obj.PlotDataViewIds)
                return;
            end
            plotIndex = get(obj.PlotViewPopup, 'Value');
            viewId = obj.PlotViewIds{plotIndex};
            if strcmp(viewId, 'all')
                return;
            end
            dataIndex = find(strcmp(viewId, obj.PlotDataViewIds), 1, 'first');
            if ~isempty(dataIndex)
                set(obj.PlotDataViewPopup, 'Value', dataIndex);
                obj.updatePlotDataTable();
            end
        end

        function updatePlotViewChoices(obj)
            views = pd_result_plot_views(obj.CurrentResult);
            labels = [{pd_ui_text('All views')}, {views.Label}];
            obj.PlotViewIds = [{'all'}, {views.Id}];
            set(obj.PlotViewPopup, 'String', labels, 'Value', 1);
            label = analysisDescription(obj.CurrentResult.analysisType);
            set(obj.PlotViewCountText, 'String', sprintf( ...
                pd_ui_text('%s / Total: %d views'), label, numel(views)));
            obj.PlotDataViewIds = {views.Id};
            set(obj.PlotDataViewPopup, 'String', {views.Label}, 'Value', 1);
            obj.PlotDataPage = 1;
            obj.PlotDataSelectedRows = [];
            obj.updatePlotDataTable();
        end

        function showPlotDataTab(obj, ~, ~)
            if isempty(obj.CurrentResult)
                obj.showError(pd_ui_text('Run an analysis before viewing plot data.'));
                return;
            end
            plotIndex = get(obj.PlotViewPopup, 'Value');
            viewId = obj.PlotViewIds{plotIndex};
            if ~strcmp(viewId, 'all')
                dataIndex = find(strcmp(viewId, obj.PlotDataViewIds), 1, 'first');
                if ~isempty(dataIndex)
                    set(obj.PlotDataViewPopup, 'Value', dataIndex);
                end
            end
            obj.updatePlotDataTable();
            set(obj.TabGroup, 'SelectedTab', obj.PlotDataTab);
        end

        function plotDataViewChanged(obj, ~, ~)
            obj.PlotDataPage = 1;
            obj.PlotDataSelectedRows = [];
            obj.updatePlotDataTable();
        end

        function updatePlotDataTable(obj)
            if isempty(obj.CurrentResult) || isempty(obj.PlotDataViewIds)
                obj.CurrentPlotData = [];
                set(obj.PlotDataTable, 'Data', cell(0, 0), 'ColumnName', {});
                set(obj.PlotDataInfoText, 'String', ...
                    pd_ui_text('Run an analysis to inspect plot data.'));
                set(obj.PlotDataPageText, 'String', pd_ui_text('Page 0 of 0'));
                set([obj.PlotDataPreviousButton,obj.PlotDataNextButton], ...
                    'Enable', 'off');
                return;
            end
            index = get(obj.PlotDataViewPopup, 'Value');
            index = max(1, min(index, numel(obj.PlotDataViewIds)));
            viewId = obj.PlotDataViewIds{index};
            obj.CurrentPlotData = pd_result_plot_data(obj.CurrentResult, viewId);
            rowCount = obj.CurrentPlotData.RowCount;
            columnCount = obj.CurrentPlotData.ColumnCount;
            pageCount = max(1, ceil(rowCount / obj.PlotDataPageSize));
            obj.PlotDataPage = max(1, min(obj.PlotDataPage, pageCount));
            firstRow = (obj.PlotDataPage - 1) * obj.PlotDataPageSize + 1;
            lastRow = min(rowCount, firstRow + obj.PlotDataPageSize - 1);
            if rowCount == 0
                pageValues = zeros(0, columnCount);
                firstRow = 0;
                lastRow = 0;
            else
                pageValues = obj.CurrentPlotData.Values(firstRow:lastRow, :);
            end
            editable = false(1, max(1, columnCount));
            if columnCount == 0
                editable = false;
            end
            set(obj.PlotDataTable, 'Data', pageValues, ...
                'ColumnName', obj.CurrentPlotData.ColumnNames, ...
                'ColumnEditable', editable);
            set(obj.PlotDataInfoText, 'String', sprintf( ...
                pd_ui_text('%s: %d rows x %d columns; showing rows %d-%d.'), ...
                obj.CurrentPlotData.Label, rowCount, columnCount, firstRow, lastRow));
            set(obj.PlotDataPageText, 'String', sprintf( ...
                pd_ui_text('Page %d of %d'), obj.PlotDataPage, pageCount));
            if obj.PlotDataPage > 1
                previousState = 'on';
            else
                previousState = 'off';
            end
            if obj.PlotDataPage < pageCount
                nextState = 'on';
            else
                nextState = 'off';
            end
            set(obj.PlotDataPreviousButton, 'Enable', previousState);
            set(obj.PlotDataNextButton, 'Enable', nextState);
            obj.PlotDataSelectedRows = [];
        end

        function changePlotDataPage(obj, delta)
            if isempty(obj.CurrentPlotData)
                return;
            end
            pageCount = max(1, ceil(obj.CurrentPlotData.RowCount / obj.PlotDataPageSize));
            obj.PlotDataPage = max(1, min(pageCount, obj.PlotDataPage + delta));
            obj.updatePlotDataTable();
        end

        function selectPlotDataRows(obj, ~, event)
            if isempty(event) || isempty(event.Indices) || isempty(obj.CurrentPlotData)
                obj.PlotDataSelectedRows = [];
                return;
            end
            pageRows = unique(event.Indices(:, 1));
            offset = (obj.PlotDataPage - 1) * obj.PlotDataPageSize;
            obj.PlotDataSelectedRows = offset + pageRows(:).';
        end

        function copyPlotData(obj, selectedOnly)
            if isempty(obj.CurrentPlotData)
                obj.showError(pd_ui_text('Run an analysis before copying plot data.'));
                return;
            end
            if selectedOnly
                rows = obj.PlotDataSelectedRows;
                if isempty(rows)
                    obj.showError(pd_ui_text('Select one or more data rows first.'));
                    return;
                end
            else
                rows = 1:obj.CurrentPlotData.RowCount;
            end
            try
                includeHeader = logical(get( ...
                    obj.PlotDataCopyHeadersCheckbox, 'Value'));
                text = pd_plot_data_text(obj.CurrentPlotData, char(9), ...
                    rows, includeHeader);
                clipboard('copy', text);
                obj.setStatus(sprintf(pd_ui_text('Copied %d plot-data rows to the clipboard.'), ...
                    numel(rows)), 'success');
            catch err
                obj.showError(err);
            end
        end

        function exportCurrentPlotData(obj, ~, ~)
            if isempty(obj.CurrentPlotData)
                obj.showError(pd_ui_text('Run an analysis before exporting plot data.'));
                return;
            end
            stepText = 'result';
            if isfield(obj.CurrentResult, 'timestep') && ...
                    isscalar(obj.CurrentResult.timestep)
                stepText = sprintf('t%g', obj.CurrentResult.timestep);
            end
            defaultName = sprintf('plotdata_%s_%s_%s.csv', ...
                obj.CurrentResult.analysisType, obj.CurrentPlotData.ViewId, stepText);
            defaultName = regexprep(defaultName, '[^A-Za-z0-9_.-]', '_');
            [name, folder] = uiputfile('*.csv', ...
                pd_ui_text('Export plotted data as CSV'), defaultName);
            if isequal(name, 0)
                return;
            end
            try
                filePath = fullfile(folder, name);
                pd_write_plot_data_csv(filePath, obj.CurrentPlotData);
                obj.setStatus([pd_ui_text('Plot data exported:'), filePath], 'success');
            catch err
                obj.showError(err);
            end
        end

        function preparePlotAxes(obj, viewIds, updateMode)
            reuse = strcmpi(updateMode, 'overlay') && ...
                numel(obj.PlotAxes) == numel(viewIds) && ...
                all(ishghandle(obj.PlotAxes));
            if reuse
                for i = 1:numel(viewIds)
                    reuse = reuse && strcmp(get(obj.PlotAxes(i), 'Tag'), ...
                        ['POST_DATA2_View_', viewIds{i}]);
                end
            end
            if reuse, return; end
            delete(get(obj.PlotPanel, 'Children'));
            count = numel(viewIds);
            [rows, columns] = pd_subplot_grid(count);
            gapX = 0.055;
            gapY = 0.075;
            marginX = 0.055;
            marginY = 0.075;
            width = (1 - 2 * marginX - (columns - 1) * gapX) / columns;
            height = (1 - 2 * marginY - (rows - 1) * gapY) / rows;
            axesList = gobjects(1, count);
            for i = 1:count
                row = floor((i - 1) / columns);
                column = mod(i - 1, columns);
                left = marginX + column * (width + gapX);
                bottom = 1 - marginY - (row + 1) * height - row * gapY;
                axesList(i) = axes('Parent', obj.PlotPanel, 'Units', 'normalized', ...
                    'Position', [left, bottom, width, height], ...
                    'Tag', ['POST_DATA2_View_', viewIds{i}]);
            end
            obj.PlotAxes = axesList;
        end

        function openAllPlotFigures(obj, ~, ~)
            if isempty(obj.CurrentResult)
                obj.showError('Run an analysis before opening result figures.');
                return;
            end
            try
                options = pd_plot_options_from_table(obj.PlotCatalog, ...
                    get(obj.PlotTable, 'Data'));
                obj.renderSeparateFigures(options, true);
                obj.setStatus(sprintf('Opened %d separate result figures.', ...
                    numel(obj.PlotFigureHandles)));
            catch err
                obj.showError(err);
            end
        end

        function renderSeparateFigures(obj, options, forceRecreate)
            views = pd_result_plot_views(obj.CurrentResult);
            viewIds = {views.Id};
            reuse = ~forceRecreate && numel(obj.PlotFigureHandles) == numel(viewIds) && ...
                numel(obj.PlotFigureAxes) == numel(viewIds) && ...
                all(ishghandle(obj.PlotFigureHandles)) && ...
                all(ishghandle(obj.PlotFigureAxes)) && ...
                isequal(obj.PlotFigureViewIds, viewIds);
            if ~reuse
                obj.closeSeparateFigures();
                count = numel(viewIds);
                obj.PlotFigureHandles = gobjects(1, count);
                obj.PlotFigureAxes = gobjects(1, count);
                for i = 1:count
                    left = 100 + 32 * mod(i - 1, 8);
                    bottom = 80 + 28 * mod(i - 1, 8);
                    obj.PlotFigureHandles(i) = figure('Name', ...
                        ['POST_DATA2 - ', views(i).Label], 'NumberTitle', 'off', ...
                        'Tag', 'POST_DATA2_ResultFigure', 'Color', 'w', ...
                        'Position', [left, bottom, 760, 560]);
                    obj.PlotFigureAxes(i) = axes('Parent', obj.PlotFigureHandles(i), ...
                        'Units', 'normalized', 'Position', [0.12, 0.12, 0.80, 0.80]);
                end
                obj.PlotFigureViewIds = viewIds;
            end
            for i = 1:numel(viewIds)
                pd_render_result(obj.PlotFigureAxes(i), obj.CurrentResult, ...
                    options, viewIds{i});
            end
        end

        function value = hasAnySeparateFigures(obj)
            value = ~isempty(obj.PlotFigureHandles) && ...
                any(ishghandle(obj.PlotFigureHandles));
        end

        function closeSeparateFigures(obj)
            if ~isempty(obj.PlotFigureHandles)
                valid = ishghandle(obj.PlotFigureHandles);
                delete(obj.PlotFigureHandles(valid));
            end
            obj.PlotFigureHandles = gobjects(1, 0);
            obj.PlotFigureAxes = gobjects(1, 0);
            obj.PlotFigureViewIds = {};
        end

        function closePlotFiguresFromButton(obj, ~, ~)
            count = nnz(ishghandle(obj.PlotFigureHandles));
            obj.closeSeparateFigures();
            obj.setStatus(sprintf(pd_ui_text('Closed %d separate result figures.'), count));
        end

        function closeApplication(obj, source, ~)
            obj.closeSeparateFigures();
            if nargin >= 2 && ishghandle(source)
                set(source, 'CloseRequestFcn', '');
                delete(source);
            end
        end

        function options = getOutputOptions(obj)
            options = pd_output_options_from_table(obj.OutputCatalog, get(obj.OutputTable, 'Data'));
        end

        function exportCurrentResult(obj, ~, ~)
            if isempty(obj.CurrentResult)
                obj.showError('Run an analysis before exporting.');
                return;
            end
            try
                obj.exportWithOptions(obj.getOutputOptions());
            catch err
                obj.showError(err);
            end
        end

        function exportWithOptions(obj, options)
            files = pd_export_result(obj.CurrentResult, obj.PlotAxes, options);
            obj.setStatus(sprintf('Exported %d file(s) to %s.', numel(files), options.Directory));
        end

        function saveConfiguration(obj, ~, ~)
            [name, folder] = uiputfile('*.mat', 'Save POST_DATA configuration', 'postdata_config.mat');
            if isequal(name, 0), return; end
            config = obj.captureConfiguration();
            save(fullfile(folder, name), 'config');
            obj.setStatus('Configuration saved.');
        end

        function config = captureConfiguration(obj)
            config = struct();
            config.schemaVersion = 2;
            config.filePath = get(obj.FileEdit, 'String');
            config.task = obj.currentAnalysisType();
            config.selectionMode = popupValue(obj.SelectionPopup);
            config.selectionValue = get(obj.SelectionEdit, 'String');
            config.slurmPath = get(obj.SlurmEdit, 'String');
            config.progressMode = popupValue(obj.ProgressPopup);
            config.resultLevel = popupValue(obj.ResultLevelPopup);
            config.computeData = get(obj.ComputeTable, 'Data');
            config.plotData = get(obj.PlotTable, 'Data');
            config.outputData = get(obj.OutputTable, 'Data');
            config.copyPlotDataHeaders = logical(get( ...
                obj.PlotDataCopyHeadersCheckbox, 'Value'));
        end

        function loadConfiguration(obj, ~, ~)
            [name, folder] = uigetfile('*.mat', 'Load POST_DATA configuration');
            if isequal(name, 0), return; end
            loaded = load(fullfile(folder, name), 'config');
            if ~isfield(loaded, 'config')
                obj.showError('Unsupported configuration file.');
                return;
            end
            try
                config = pd_upgrade_app_config(loaded.config);
            catch err
                obj.showError(err);
                return;
            end
            obj.applyConfiguration(config);
            obj.setStatus('Configuration loaded.');
        end

        function applyConfiguration(obj, config)
            set(obj.FileEdit, 'String', config.filePath);
            obj.setAnalysisType(config.task);
            obj.taskChanged([], []);
            setPopupValue(obj.SelectionPopup, config.selectionMode);
            set(obj.SelectionEdit, 'String', config.selectionValue);
            set(obj.SlurmEdit, 'String', config.slurmPath);
            setPopupValue(obj.ProgressPopup, config.progressMode);
            setPopupValue(obj.ResultLevelPopup, config.resultLevel);
            if size(config.computeData, 1) == numel(obj.OptionCatalog)
                computeData = pd_normalize_catalog_table_data( ...
                    obj.OptionCatalog, config.computeData);
                set(obj.ComputeTable, 'Data', computeData);
                obj.updateDependencyDisplay();
            end
            if size(config.plotData, 1) == numel(obj.PlotCatalog)
                plotData = pd_normalize_catalog_table_data( ...
                    obj.PlotCatalog, config.plotData);
                set(obj.PlotTable, 'Data', plotData);
            end
            if size(config.outputData, 1) == numel(obj.OutputCatalog)
                outputData = pd_normalize_catalog_table_data( ...
                    obj.OutputCatalog, config.outputData);
                set(obj.OutputTable, 'Data', outputData);
            end
            set(obj.PlotDataCopyHeadersCheckbox, 'Value', ...
                double(logical(config.copyPlotDataHeaders)));
            obj.selectionModeChanged([], []);
            obj.updateRunButtonForParameters();
        end

        function setStatus(obj, message, state)
            if nargin < 3, state = 'normal'; end
            if ishghandle(obj.StatusText)
                switch lower(state)
                    case 'busy'
                        color = [1.0, 0.96, 0.78];
                    case 'success'
                        color = [0.86, 0.97, 0.86];
                    case 'error'
                        color = [1.0, 0.84, 0.84];
                    otherwise
                        color = 'white';
                end
                set(obj.StatusText, 'String', message, 'BackgroundColor', color);
                drawnow;
            end
        end

        function showError(obj, errorValue)
            if isa(errorValue, 'MException')
                message = localizeValidationMessage(errorValue.message);
                identifier = errorValue.identifier;
                diagnostic = getReport(errorValue, 'extended', ...
                    'hyperlinks', 'off');
                logFile = pd_default_log_file();
                pd_log(logFile, 'error', 'app', diagnostic);
                if isempty(identifier), identifier = 'postdata:UnknownError'; end
                dialogMessage = {message, [pd_ui_text('Error identifier:'), identifier], ...
                    [pd_ui_text('Detailed log:'), logFile]};
            else
                message = pd_to_char(errorValue);
                dialogMessage = {message};
            end
            obj.setStatus([pd_ui_text('Error:'), message], 'error');
            errordlg(dialogMessage, 'POST_DATA2');
        end
    end
end

function tableHandle = makeOptionTable(parent)
    tableHandle = uitable('Parent', parent, 'Units', 'normalized', ...
        'Position', [0.02, 0.11, 0.96, 0.86], ...
        'ColumnName', {pd_ui_text('Parameter'),pd_ui_text('Value'),pd_ui_text('Description')}, ...
        'ColumnEditable', [false true false], ...
        'ColumnWidth', {190 180 430}, 'Data', cell(0, 3));
end

function addLabel(parent, label, y)
    uicontrol('Parent', parent, 'Style', 'text', 'String', label, ...
        'Units', 'normalized', 'Position', [0.04, y, 0.92, 0.035], ...
        'HorizontalAlignment', 'left');
end

function control = addEdit(parent, value, position)
    control = uicontrol('Parent', parent, 'Style', 'edit', 'String', value, ...
        'Units', 'normalized', 'Position', position, ...
        'HorizontalAlignment', 'left', 'BackgroundColor', 'white');
end

function control = addPopup(parent, values, position, callback)
    control = uicontrol('Parent', parent, 'Style', 'popupmenu', 'String', values, ...
        'Units', 'normalized', 'Position', position, 'BackgroundColor', 'white');
    if ~isempty(callback), set(control, 'Callback', callback); end
end

function control = addButton(parent, label, position, callback)
    control = uicontrol('Parent', parent, 'Style', 'pushbutton', 'String', label, ...
        'Units', 'normalized', 'Position', position, 'Callback', callback);
end

function value = popupValue(control)
    values = get(control, 'String');
    index = get(control, 'Value');
    if ischar(values)
        value = strtrim(values(index, :));
    else
        value = values{index};
    end
end

function setPopupValue(control, value)
    values = get(control, 'String');
    if ischar(values), values = cellstr(values); end
    index = find(strcmpi(value, values), 1, 'first');
    if isempty(index)
        error('postdata:BadPopupValue', 'Unsupported value: %s', value);
    end
    set(control, 'Value', index);
end

function data = replaceTableOption(data, name, value)
    row = find(strcmp(name, data(:, 1)), 1, 'first');
    if isempty(row)
        error('postdata:MissingPlotPresetOption', ...
            'Plot preset refers to unknown option: %s.', name);
    end
    data{row, 2} = value;
end

function text = analysisDescription(task)
    switch lower(strtrim(pd_to_char(task)))
        case 'chunk'
            text = pd_ui_text('1D/2D field analysis');
        case 'network2d'
            text = pd_ui_text('2D pore network analysis');
        case 'cluster'
            text = pd_ui_text('Particle/cluster analysis');
        case 'vx'
            text = pd_ui_text('mass-v analysis');
        case 'massx'
            text = pd_ui_text('SPH mass-x (m-x) analysis');
        otherwise
            text = pd_to_char(task);
    end
end

function text = computeHeaderText(task, count)
    if strcmp(task, 'massx')
        text = sprintf(pd_ui_text(['SPH mass-x (m-x): InitialDensity and ', ...
            'ParticleSpacing are required; bin1d also requires ', ...
            'TransverseWidth. (%d parameters)']), count);
    else
        text = sprintf(pd_ui_text('%s: %d editable calculation parameters'), ...
            analysisDescription(task), count);
    end
end

function message = localizeValidationMessage(message)
    message = strrep(pd_to_char(message), 'A value is required.', ...
        pd_ui_text('A value is required.'));
end

function position = defaultWindowPosition()
    screen = get(0, 'ScreenSize');
    width = max(520, min(1320, screen(3) - 40));
    height = max(440, min(780, screen(4) - 80));
    left = max(1, screen(1) + (screen(3) - width) / 2);
    bottom = max(1, screen(2) + (screen(4) - height) / 2);
    position = round([left, bottom, width, height]);
end

function initial = browseStart(currentValue, pattern)
    currentValue = strtrim(pd_to_char(currentValue));
    if exist(currentValue, 'file')
        folder = fileparts(currentValue);
    elseif exist(currentValue, 'dir')
        folder = currentValue;
    else
        folder = pwd;
    end
    initial = fullfile(folder, pattern);
end
