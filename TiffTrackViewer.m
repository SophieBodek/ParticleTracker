classdef TiffTrackViewer < handle
    %Class to view tiffs using memory mapping
    %tv=TiffTrackViewer(input,n_ch,vtracks,min_track_length);
    %
    %input = .tif file path, .bin file path, folder (for multi-file tif stack)
    %or matrix
    %
    %n_ch = number of channels in the file (default behavior: will try to
    %determine from the file, if possible; if not, 1).
    %
    %vtracks = struct array with particle tracking data (optional)
    %
    %min_track_length = minimum track length to display (optional, default = 1)
    %
    %GUI will popup and has intuitive controls.
    %
    %Controls:
    %
    %Play / Pause button - plays videos starting at current frame
    %scroll - scroll through video frames
    %Current Frame input - type in frame that you want to navigate to
    %Set FPS - set video speed in frames per second
    %Tracks - toggle particle tracks on and off; can click on tracks to display track info
    %Max P - Maximum projection with a little noise removal first
    %Mean P - Mean projection (i.e., background image)
    %
    %Projections are currently capped at taking 30 seconds max to avoid
    %stalling Matlab.
    %
    %This code is adapted from TiffViewer/Movieviewer, written by
    %Joseph M. Stujenske (under GPL-3.0 copyright license) and is 
    %avaible in the Matlab_FastTiffReadWrite Github repository 
    %(accessed 15 July 2025). This current version was written with
    %the assistance of generative AI tools.
    %
    properties
        filename                % Name of file
        figure
        n_ch = 1;               % Number of channels (defaults to 1)
        numFrames               % Number of frames
        memmap_data
        memmap_matrix_data
        fps                     % Frames per second
        map_type
        max_time = 20;
        type
        width
        height
        CurrFrame = 1;
        vtracks = [];           % Particle tracking data
        show_tracks = false;    % Toggle for showing tracks
        track_handles = {};     % Store line handles for tracks
        future_track_handles = {}; % Store handles for future trajectory lines
        min_track_length = 1;   % Minimum track length filter
        filtered_vtracks = [];  % Filtered track data
    end
    
    properties(Hidden = true)
        ax
        memmap
        memmap_matrix
        listener
        image_object
        dim_text                % Text object for displaying dimensions
    end
    
    methods
        function obj = TiffTrackViewer(filename, n_ch, vtracks, min_track_length)
            if nargin < 4, obj.min_track_length = 1; else obj.min_track_length = min_track_length; end % NEW
            if nargin < 3, obj.vtracks = []; else obj.vtracks = vtracks; end % NEW
            if nargin < 2, obj.n_ch = []; else obj.n_ch=n_ch;end
            if nargin < 1 || isempty(filename)
                filename = uigetfile({'*.tif;*.tiff'});
                if isempty(obj.n_ch)
                    userinputs = inputdlg({'Number of channels:'});
                    obj.n_ch = str2double(userinputs{1});
                end
            end

            % Initialize future track handles
            obj.future_track_handles = {};

            % Filter tracks by length
            obj.filterTracksByLength();
            
            if ischar(filename)
                obj.resolveFileName(filename);
                [~,~,ext]=fileparts(filename);
                if ~isempty(ext)
                info = obj.readTiffInfo;
                else
                    info=[];
                end
                obj.setNumFrames(info);
                uneven_flag = obj.checkUnevenFlag(info);
                obj.setupMemoryMapping(filename, obj.n_ch, info, uneven_flag);
                if isempty(info)
                    if iscell(obj.filename)
                    info = obj.readTiffInfo(obj.filename{1});
                    end
                end
                dims = obj.getDimensions(info);
                file=filename;
            else
                if isempty(obj.n_ch), obj.n_ch = 1; end
                obj.setupFromMatrix(filename, obj.n_ch);
                file='Matrix';
                dims=size(filename);
                ext=[];
            end
            obj.setupFigure(file, ext);
            obj.setupAxes(dims, obj.n_ch);
            obj.setupSlider();
            obj.setupTimer();
            obj.addDimensionText(dims); % Add dimension text display
            obj.displayFrame(1);
        end
        
        function displayFrame(obj,opt)
            if nargin < 2, opt = 0; end
            obj.validateFrame();
            for a = 1:obj.n_ch
                if obj.map_type == "mem"
                    obj.updateFrameFromMemory(a, opt);
                else
                    obj.updateFrameFromFile(a, opt);
                end
            end
            obj.updateTrackDisplay(); % Update track overlay
        end

        function allimages = mm_proj(obj, type)
            allimages = mm_proj([], [], obj, type);
        end

    % Toggle track display
    function toggleTracks(obj)
        % Check if track data is available
        if isempty(obj.vtracks)
            fprintf('No track data provided\n');
            msgbox('No track data provided. Please load vtracks when calling the viewer.', ...
                   'Track Display', 'warn');
            return;
        end
        
        % Check if filtered tracks are available after length filtering
        if isempty(obj.filtered_vtracks)
            fprintf('No tracks meet the minimum length criteria\n');
            fprintf('Original tracks: %d, Minimum length: %d, Tracks after filtering: %d\n', ...
                    length(obj.vtracks), obj.min_track_length, length(obj.filtered_vtracks));
            
            msg = sprintf(['No tracks meet the minimum length criteria.\n\n' ...
                          'Original tracks: %d\n' ...
                          'Minimum length: %d frames\n' ...
                          'Tracks after filtering: %d\n\n' ...
                          'Try reducing the minimum track length.'], ...
                          length(obj.vtracks), obj.min_track_length, length(obj.filtered_vtracks));
            
            msgbox(msg, 'Track Display', 'warn');
            return;
        end
        
        % Toggle track display
        obj.show_tracks = ~obj.show_tracks;
        if ~obj.show_tracks
            obj.clearTracks();
            fprintf('Track display turned OFF\n');
        else
            fprintf('Track display turned ON\n');
            fprintf('Displaying %d tracks (filtered from %d original tracks)\n', ...
                    length(obj.filtered_vtracks), length(obj.vtracks));
        end
        obj.updateTrackDisplay();
    end

    % Clear all track lines from display
    function clearTracks(obj)
        for ch = 1:obj.n_ch
            if ~isempty(obj.track_handles) && length(obj.track_handles) >= ch
                delete(obj.track_handles{ch}(ishandle(obj.track_handles{ch})));
                obj.track_handles{ch} = [];
            end
        end
    end

    % Filter tracks by minimum length
    function filterTracksByLength(obj)
        if isempty(obj.vtracks)
            obj.filtered_vtracks = [];
            return;
        end
        
        % Calculate track lengths and filter
        track_lengths = arrayfun(@(x) x.len, obj.vtracks);
        valid_tracks = track_lengths >= obj.min_track_length;
        
        obj.filtered_vtracks = obj.vtracks(valid_tracks);
        
        % Display filtering results
        fprintf('Track filtering: %d/%d tracks retained (min length: %d)\n', ...
            sum(valid_tracks), length(obj.vtracks), obj.min_track_length);
    end
    
    % Update minimum track length filter
    function setMinTrackLength(obj, min_length)
        obj.min_track_length = min_length;
        obj.filterTracksByLength();
        if obj.show_tracks
            obj.updateTrackDisplay();
        end
    end

    function updateTrackDisplay(obj)
        if isempty(obj.filtered_vtracks) || ~obj.show_tracks % Use filtered_vtracks instead of vtracks
            obj.clearFutureTrajectories(); % Clear future trajectories when tracks are hidden
            return;
        end
        
        obj.clearTracks();
        obj.clearFutureTrajectories(); % Clear future trajectories when updating normal tracks
        
        colorlist = lines(7);
        Ncolors = size(colorlist, 1);
        ntracks = numel(obj.filtered_vtracks); % Use filtered_vtracks
        
        % Calculate frame range for each track
        framerange = NaN(ntracks, 2);
        for jj = 1:ntracks
            framerange(jj, 1) = obj.filtered_vtracks(jj).T(1); % Use filtered_vtracks
            framerange(jj, 2) = obj.filtered_vtracks(jj).T(end); % Use filtered_vtracks
        end
        
        % Find tracks that should be visible in current frame
        ind = find((framerange(:, 1) <= obj.CurrFrame) & (framerange(:, 2) >= obj.CurrFrame));
        
        % Plot tracks for each channel
        for ch = 1:obj.n_ch
            obj.track_handles{ch} = [];
            
            for jj = 1:numel(ind)
                track_idx = ind(jj);
                col = colorlist(mod(track_idx - 1, Ncolors) + 1, :);
                
                % Plot from beginning of track to current frame
                indt = 1:(obj.CurrFrame - framerange(track_idx, 1) + 1);
                
                if ~isempty(indt) && indt(end) <= length(obj.filtered_vtracks(track_idx).X)
                    hold(obj.ax{ch}, 'on');
                    h = plot(obj.ax{ch}, ...
                        obj.filtered_vtracks(track_idx).X(indt), ...
                        obj.filtered_vtracks(track_idx).Y(indt), ...
                        '-', 'Color', col, 'LineWidth', 1);
                    

                    % Determine track ID: use ID field if it exists, otherwise use index
                    if isfield(obj.filtered_vtracks, 'ID') && ~isempty(obj.filtered_vtracks(track_idx).ID)
                        track_id = obj.filtered_vtracks(track_idx).ID;
                    else
                        track_id = track_idx; % Use index (of filtered vtracks) as fallback
                    end
                    
                    % Store track information in the line's UserData
                    set(h, 'UserData', struct('TrackID', track_id, 'TrackIndex', track_idx));
                    
                    % Add click callback
                    set(h, 'ButtonDownFcn', @(src, event) obj.onTrackClick(src, event));
                    
                    % Make the line more clickable by increasing hit test area
                    set(h, 'LineWidth', 2, 'PickableParts', 'visible');
                    
                    obj.track_handles{ch} = [obj.track_handles{ch}, h];
                    hold(obj.ax{ch}, 'off');
                end
            end
        end
    end

    % NEW METHOD: Show future trajectory as dashed line
    function showFutureTrajectory(obj, track_info, track_idx)
        % Clear any existing future trajectories
        obj.clearFutureTrajectories();

        % Initialize future_track_handles if needed
        if isempty(obj.future_track_handles)
            obj.future_track_handles = cell(1, obj.n_ch);
        end
        
        % Ensure we have enough cells for all channels
        while length(obj.future_track_handles) < obj.n_ch
            obj.future_track_handles{end+1} = [];
        end
            
        % Get track timing information
        track_frames = track_info.T;
        current_frame_idx = find(track_frames == obj.CurrFrame);
        
        % If current frame is not exactly on a track frame, find the next frame
        if isempty(current_frame_idx)
            future_frame_idx = find(track_frames > obj.CurrFrame, 1, 'first');
        else
            future_frame_idx = current_frame_idx + 1;
        end
        
        % Only show future if there are future points
        if ~isempty(future_frame_idx) && future_frame_idx <= length(track_info.X)
            future_indices = future_frame_idx:length(track_info.X);
            future_X = track_info.X(future_indices);
            future_Y = track_info.Y(future_indices);
            
            % Get color for this track (same as the solid line)
            colorlist = lines(7);
            Ncolors = size(colorlist, 1);
            col = colorlist(mod(track_idx - 1, Ncolors) + 1, :);
            
            % Plot future trajectory on each channel
            for ch = 1:obj.n_ch
                if ~isempty(future_X) && ~isempty(future_Y)
                    hold(obj.ax{ch}, 'on');
                    h_future = plot(obj.ax{ch}, future_X, future_Y, ...
                                   '--', 'Color', col, 'LineWidth', 2, ...
                                   'DisplayName', 'Future Trajectory');
                    obj.future_track_handles{ch} = [obj.future_track_handles{ch}, h_future];
                    hold(obj.ax{ch}, 'off');

                    % Initialize channel cell if needed
                    if isempty(obj.future_track_handles{ch})
                        obj.future_track_handles{ch} = [];
                    end
                end
                obj.future_track_handles{ch} = [obj.future_track_handles{ch}, h_future];
                hold(obj.ax{ch}, 'off');
            end
            
            fprintf('Future trajectory displayed (%d future points)\n', length(future_indices));
        else
            fprintf('No future trajectory to display (track ends at or before current frame)\n');
        end
    end

    % NEW METHOD: Clear future trajectory lines
    function clearFutureTrajectories(obj)
        if isempty(obj.future_track_handles)
            return;
        end
        
        for ch = 1:min(length(obj.future_track_handles), obj.n_ch)
            if ~isempty(obj.future_track_handles{ch})
                % Delete valid handles
                valid_handles = obj.future_track_handles{ch}(ishandle(obj.future_track_handles{ch}));
                if ~isempty(valid_handles)
                    delete(valid_handles);
                end
                obj.future_track_handles{ch} = [];
            end
        end
end

    % Handle track click events
    function onTrackClick(obj, src, event)
        % Reset all track line widths to normal and clear any future trajectories
        obj.clearFutureTrajectories();
        for ch = 1:obj.n_ch
            if ~isempty(obj.track_handles{ch})
                set(obj.track_handles{ch}, 'LineWidth', 2);
            end
        end
        
        % Highlight the clicked track
        set(src, 'LineWidth', 4);
        
        % Get track information from the clicked line
        trackData = get(src, 'UserData');
        track_id = trackData.TrackID;
        track_idx = trackData.TrackIndex;
        
        % Get additional track information
        track_info = obj.filtered_vtracks(track_idx);
        track_length = track_info.len;
        start_frame = track_info.T(1);
        end_frame = track_info.T(end);
        
        % Determine if we're using actual ID or index
        if isfield(obj.filtered_vtracks, 'ID') && ~isempty(obj.filtered_vtracks(track_idx).ID)
            id_type = 'ID (Original)';
        else
            id_type = 'Index (Filtered)';
        end

        % Print to command window
        fprintf('Clicked Track %s: %d (Length: %d frames, Frames %d-%d)\n', ...
                id_type, track_id, track_length, start_frame, end_frame);

        % Show future trajectory if current frame is before track end
        if obj.CurrFrame < end_frame
            obj.showFutureTrajectory(track_info, track_idx);
        end
        
        % Create custom dialog with plot option
        obj.showTrackInfoDialog(track_info, track_id, track_idx, id_type, track_length, start_frame, end_frame);
    end

    % NEW METHOD: Show track information dialog with plot option
    function showTrackInfoDialog(obj, track_info, track_id, track_idx, id_type, track_length, start_frame, end_frame)
        % Calculate quick statistics for display
        % if isfield(track_info, 'U') && isfield(track_info, 'V')
        %     U = track_info.U;
        %     V = track_info.V;
        %     vel = sqrt(U.^2 + V.^2);
        %     mean_speed = nanmean(vel);
        %     max_speed = nanmax(vel);
        %     speed_info = sprintf('\nMean Speed: %.2f px/frame\nMax Speed: %.2f px/frame', mean_speed, max_speed);
        % else
        %     speed_info = '\nVelocity data not available';
        % end
        
        % Create the dialog message
        msg = sprintf(['Track %s: %d\n' ...
                      'Track Length: %d frames\n' ...
                      'Start Frame: %d\n' ...
                      'End Frame: %d\n' ...
                      'Current Frame: %d%s\n\n' ...
                      'Would you like to plot track\n' ...
                      'position and velocity timeseries?'], ...
                      id_type, track_id, track_length, start_frame, end_frame, obj.CurrFrame);%, speed_info);
        
        % Create custom dialog with Yes/No buttons
        choice = questdlg(msg, ...
                         'Track Information', ...
                         'Plot Trajectories', 'Close', 'Close');
        
        % Handle the user's choice
        switch choice
            case 'Plot Trajectories'
                obj.plotTrackDetails(track_info, track_id, id_type);
            case 'Close'
                % Do nothing, just close the dialog
                return;
        end
    end
    
    % NEW METHOD: Create detailed trajectory and velocity plots
    function plotTrackDetails(obj, track_info, track_id, id_type)
        % Extract track data
        len = track_info.len;
        X = track_info.X;
        Y = track_info.Y;
        T = track_info.T;
        
        % Calculate displacements
        dX = [nan diff(X)];
        dY = [nan diff(Y)];
        displacement = sqrt(dX.^2 + dY.^2);
        
        % Extract or calculate velocities
        if isfield(track_info, 'U') && isfield(track_info, 'V')
            U = track_info.U;
            V = track_info.V;
        else
            % Calculate velocities if not available (backward difference)
            dt = 1; %timestep = 1 frame
            U = [nan diff(X)./dt]; % X velocity - pad with NaN at start
            V = [nan diff(Y)./dt]; % Y velocity - pad with NaN at start
            fprintf('Note: Velocities calculated from positions (U and V fields not found)\n');
        end
        
        vel = sqrt(U.^2 + V.^2); % speed
        
        % Calculate statistics
        mean_speed = nanmean(vel);
        max_speed = nanmax(vel);
        total_distance = nansum(displacement);
        net_displacement = sqrt((X(end)-X(1))^2 + (Y(end)-Y(1))^2);
        
        % Print statistics
        fprintf('Track %d Statistics:\n', track_id);
        fprintf('Mean speed: %.2f px/frame\n', mean_speed);
        fprintf('Max speed: %.2f px/frame\n', max_speed);
        fprintf('Total distance: %.2f px\n', total_distance);
        fprintf('Net displacement: %.2f px\n', net_displacement);
        
        % Create the detailed plot
        figure('Name', sprintf('Track %s %d Details', id_type, track_id), ...
               'NumberTitle', 'off');
        tiledlayout(5,1, 'TileSpacing', 'compact', 'Padding', 'compact');
        set(gcf, 'Position', get(0, 'Screensize'));
        
        % X position and displacement
        nexttile
        yyaxis left
        plot(T, X, 'k-', 'LineWidth', 2, 'DisplayName', 'X Position [px]');
        hold on
        ylabel('X Position [px]', 'FontSize', 12);
        set(gca, 'YColor', 'k');
        yyaxis right
        plot(T, dX, 'Color', '#1ad1ff', 'LineWidth', 1.5, 'DisplayName', 'X Displacement [px]');
        ylabel('X Displacement [px]', 'FontSize', 12);
        set(gca, 'YColor', '#1ad1ff');
        title(sprintf('Track %s %d: X Position & Displacement', id_type, track_id), 'FontSize', 12);
        grid on; grid minor;
        
        % Velocity, X direction
        nexttile
        plot(T, U, 'Color', '#e60000', 'LineWidth', 1.5, 'DisplayName', 'Velocity, U');
        ylabel('Velocity [px/frame]', 'FontSize', 12);
        set(gca, 'YColor', '#e60000');
        title('X Velocity', 'FontSize', 12);
        yline(nanmean(U), '--k', sprintf('Mean [px/frame]: %.2f', nanmean(U)));
        grid on; grid minor;
        
        % Y position and displacement
        nexttile
        yyaxis left
        plot(T, Y, 'k-', 'LineWidth', 2, 'DisplayName', 'Y Position [px]');
        hold on
        ylabel('Y Position [px]', 'FontSize', 12);
        set(gca, 'YColor', 'k');
        yyaxis right
        plot(T, dY, 'Color', '#1ad1ff', 'LineWidth', 1.5, 'DisplayName', 'Y Displacement [px]');
        ylabel('Y Displacement [px]', 'FontSize', 12);
        set(gca, 'YColor', '#1ad1ff');
        title('Y Position & Displacement', 'FontSize', 12);
        grid on; grid minor;
        
        % Velocity, Y direction
        nexttile
        plot(T, V, 'Color', '#e60000', 'LineWidth', 1.5, 'DisplayName', 'Velocity, V');
        ylabel('Velocity [px/frame]', 'FontSize', 12);
        set(gca, 'YColor', '#e60000');
        title('Y Velocity', 'FontSize', 12);
        yline(nanmean(V), '--k', sprintf('Mean [px/frame]: %.2f', nanmean(V)));
        grid on; grid minor;
        
        % Net displacement and speed
        nexttile
        yyaxis left
        plot(T, displacement, 'Color', '#0000cc', 'LineWidth', 1.5, 'DisplayName', 'Net Displacement [px]');
        hold on
        ylabel('Net Displacement [px]', 'FontSize', 12);
        set(gca, 'YColor', '#0000cc');
        yyaxis right
        plot(T, vel, 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'Speed [px]');
        ylabel('Speed [px/frame]', 'FontSize', 12);
        set(gca, 'YColor', 'r');
        yline(mean_speed, '--r', sprintf('Mean [px/frame]: %.2f', mean_speed));
        title('Net Displacement & Speed', 'FontSize', 12);
        grid on; grid minor;
        
        xlabel('Time [frame no.]', 'FontSize', 12);
        
        % % Add text box with statistics
        % annotation('textbox', [0.02, 0.02, 0.3, 0.15], ...
        %            'String', sprintf(['Track %s: %d\n' ...
        %                             'Length: %d frames\n' ...
        %                             'Mean Speed: %.2f px/frame\n' ...
        %                             'Max Speed: %.2f px/frame\n' ...
        %                             'Total Distance: %.2f px\n' ...
        %                             'Net Displacement: %.2f px'], ...
        %                             id_type, track_id, len, mean_speed, max_speed, ...
        %                             total_distance, net_displacement), ...
        %            'FontSize', 10, ...
        %            'BackgroundColor', 'white', ...
        %            'EdgeColor', 'black');
    end

end
    
    methods (Access = private)
        function resolveFileName(obj, filename)
            [folder, file, ext] = fileparts(filename);
            if isempty(ext)
                temp = dir(fullfile(folder, file, '*.tif*'));
                obj.filename = arrayfun(@(x) fullfile(x.folder, x.name), temp, 'UniformOutput', false);
            elseif any(strcmp(ext, {'.tif', '.tiff', '.bin'}))
                obj.filename = filename;
            else
                obj.filename = [];
            end
        end

        function setNumFrames(obj,info)
            if isempty(info)
                obj.numFrames=length(obj.filename);
            else
                obj.numFrames=length(info);
            end
        end
        
        function info = readTiffInfo(obj,filename)
            if nargin<2
                filename=obj.filename;
            end
            [folder,file,ext]=fileparts(filename);
            if any(strcmp(ext, {'.tif', '.tiff'}))
                try
                info = readtifftags(filename);
                catch
                    info=imfinfo(filename);
                    info.ImageWidth=info.Width;
                    info.ImageHeight=info.Height;
                end
            else
                info = [];
            end
        end
        
        function uneven_flag = checkUnevenFlag(~, info)
            if length(info) > 3 && info(2).StripOffsets(1) - info(1).StripOffsets(1) ~= info(3).StripOffsets(1) - info(2).StripOffsets(1)
                uneven_flag = 1;
            elseif isempty(info) || info(1).Compression~=1
                uneven_flag = 1;
            else
                uneven_flag = 0;
            end
        end
        
        function setupMemoryMapping(obj, filename, n_ch, info, uneven_flag)
            [~,~,ext]=fileparts(filename);
            if strcmp(ext, '.bin') || ~uneven_flag
                obj.map_type = 'mem';
                obj.handleMemoryMapping(filename, n_ch, info);
                    if exist('info','var') && ~isempty(info)
                        if isfield(info,'GapBetweenImages') && info(1).GapBetweenImages==0
                            obj.n_ch=length(fieldnames(obj.memmap_data));
                        elseif length(info)>2
                            obj.n_ch=length(fieldnames(obj.memmap_data))/2;
                        else
                            obj.n_ch=length(fieldnames(obj.memmap_data));
                        end
                    else
                        if isempty(n_ch)
                            n_ch=1;
                        end
                        obj.n_ch=n_ch;
                    end
            else
                obj.map_type = 'file';
                obj.setupFileMapping(filename, info, n_ch);
                if ~isempty(n_ch)
                    obj.n_ch=n_ch;
                    else
                        n_ch=1;
                        obj.n_ch=1;
                end
            end
            n_ch=obj.n_ch;
        end

        function handleMemoryMapping(obj, filename, n_ch, info)
            [~,~,ext]=fileparts(filename);
            switch ext
                case {'.tif', '.tiff', []}
                    obj.setupTiffMapping(filename, n_ch, info);
                    if ~iscell(obj.memmap_data)
                    obj.memmap_data=obj.memmap.Data;
                    obj.map_type = 'file';
                    end
                    obj.type='tif';
                case '.bin'
                    obj.setupBinaryMapping(filename, n_ch);
                    obj.type='binary';
                otherwise
                    obj.map_type = 'file';
            end
        end
        
        function setupFileMapping(obj, filename, info, n_ch)
            % obj.memmap = set_up_file(filename, info, n_ch);
            % obj.memmap_data = obj.memmap.Data;
        end
        
        function setupFromMatrix(obj, filename, n_ch)
            % filename = permute(filename, [2 1 3]);
            [height, width, obj.numFrames] = size(filename);
            % datavals = obj.processChannels(filename, height, width, n_ch);
            % obj.memmap_data = cell2struct(datavals, obj.generateChannelNames(n_ch), 2)';
            obj.memmap_matrix_data = filename;
            obj.map_type = 'mem';
            obj.type='matrix';
        end

        % function datavals = processChannels(~, filename, height, width, n_ch)
        %     datavals = [];
        %     for ch_rep = 1:n_ch
        %         datavals = cat(2, datavals, mat2cell(filename(:, :, ch_rep:n_ch:end), height, width, ones(size(filename, 3) / n_ch, 1)));
        %     end
        % end
        
        function ch_names = generateChannelNames(~, n_ch)
            ch_names = arrayfun(@(x) ['channel', num2str(x)], 1:n_ch, 'UniformOutput', false);
        end
        
        function setupFigure(obj, file, ext)
            obj.figure = figure('Units', 'normalized', 'Position', [0 0 1 1], ...
                'AutoResizeChildren', 'off', ...
                'CloseRequestFcn', @(x, event) obj.closeFigure(x, event), ...
                'Name', [file, ext], 'NumberTitle', 'off', 'MenuBar', 'none', ...
                'ToolBar', 'none');
        end
        
        function setupAxes(obj, dims, n_ch)
            obj.n_ch = n_ch;
            height=dims(1);
            width=dims(2);
            obj.height = height;  % Store dimensions as object properties
            obj.width = width;    % Store dimensions as object properties
            for rep = 1:obj.n_ch
                obj.ax{rep} = axes('Units', 'normalized', 'Parent', obj.figure, 'Position', [0+(rep-1)*.5 0 .5 .89], 'XTick', [], 'YTick', []);
                set(obj.ax{rep}, 'XLim', [0.5 width+0.5], 'YLim', [0.5 height+0.5]);%'XLim', [1 width], 'YLim', [1 height]);
                axis(obj.ax{rep}, 'image'); % This maintains the correct aspect ratio
            end
            linkaxes(cat(1, obj.ax{:}));
            % Initialize future track handles for the number of channels
            obj.future_track_handles = cell(1, obj.n_ch); 
        end

        function dims = getDimensions(obj, info)
            if ~isempty(info)
            dims = [info(1).ImageHeight info(1).ImageWidth];
            else
                dims=[obj.height obj.width];
            end
        end

        % Add text displaying image dimensions in bottom right corner
        function addDimensionText(obj, dims)
            % Add text displaying image dimensions in bottom right corner
            if length(dims) >= 2
                height = dims(1);
                width = dims(2);
                dim_string = sprintf('%d x %d x %d', width, height, obj.numFrames);
            else
                dim_string = 'Dimensions unknown';
            end
            
            % Create text object in the bottom right corner of the figure
            obj.dim_text = uicontrol('Style', 'text', ...
                'String', dim_string, ...
                'Units', 'normalized', ...
                'Position', [0.75, 0.01, 0.24, 0.05], ...
                'BackgroundColor', [0.9 0.9 0.9], ...
                'FontSize', 10, ...
                'HorizontalAlignment', 'right', ...
                'Parent', obj.figure);
        end

        function setupSlider(tv)
        data=guidata(tv.figure);
        data.h.slide = uicontrol('style','slider','units','normalized','position',[0.05 .92 .5 .05],'Parent',tv.figure,'Max',tv.numFrames,'Min',1,'Value',1,'SliderStep',[1, 1] / (max(tv.numFrames,2) - 1));
        data.h.edit = uicontrol('style','edit','units','normalized','position',[.57 .92 .05 .05],'Parent',tv.figure,'Max',1,'Min',1,'String',num2str(1),'callback',{@(hObject, event) makeplot2(hObject, event,tv)});
        data.h.play = uicontrol('style','pushbutton','units','normalized','position',[0 .92 .05 .05],'String','>','callback',{@(hObject,event) play_but_down(hObject,event,tv)});
        data.h.setfps = uicontrol('style','pushbutton','units','normalized','position',[.65 .92 .1 .05],'String','Set FPS','callback',{@(hObject,event) fps_but_down(hObject,event,tv)});
        data.h.tracks = uicontrol('style','pushbutton','units','normalized','position',[.75 .92 .05 .05],'String','Tracks','callback',{@(hObject,event) tv.toggleTracks()}); % Track toggle button
        data.h.maxp = uicontrol('style','pushbutton','units','normalized','position',[.8 .92 .1 .05],'String','Max P','callback',{@(hObject,event) mm_proj(hObject,event,tv,'max')});
        data.h.meanp = uicontrol('style','pushbutton','units','normalized','position',[.9 .92 .1 .05],'String','Mean P','callback',{@(hObject,event) mm_proj(hObject,event,tv,'mean')});

        if tv.numFrames==1
            data.h.slide.Visible=false;
            data.h.play.Visible=false;
            data.h.setfps.Visible=false;
        end
        guidata(tv.figure,data);
        tv.listener=addlistener(data.h.slide,'ContinuousValueChange',@(hObject, event) makeplot(hObject, event,tv));
        end

        function setupTimer(obj)
            data = guidata(obj.figure);
            data.increment = 1;
            obj.CurrFrame = 1;
            obj.fps = 30;
            data.timer = timer('ExecutionMode', 'fixedRate', 'TimerFcn', {@play_vid, obj}, 'Period', max(round(1 / obj.fps, 3), .001));
            guidata(obj.figure, data);
        end

        function closeFigure(obj, ~, ~)
            data = guidata(obj.figure);
            try;delete(obj.listener);end
            try;stop(data.timer);end
            try;obj.cleanupMemoryMapping();end
            delete(obj.figure);
        end
        
        function cleanupMemoryMapping(obj)
            obj.memmap_data = [];
            obj.memmap_matrix_data = [];
            obj.memmap = [];
            obj.memmap_matrix = [];
        end
        
        function validateFrame(obj)
            if nargin < 2 || isempty(frame)
                obj.CurrFrame = min(max(obj.CurrFrame, 1), obj.numFrames);
            end
            if obj.CurrFrame > obj.numFrames
                obj.CurrFrame = obj.numFrames;
            end
        end

        function updateFrameFromMemory(obj, channel, opt)
            data=guidata(obj.figure);
            dataField = ['channel', num2str(channel)];
            frame=obj.CurrFrame;
            
            if isempty(obj.memmap_matrix_data)
            imgData = obj.memmap_data(frame).(dataField);
            else
            [y_len, x_len] = size(obj.memmap_matrix_data, 1:2);
            if strcmp(obj.type,'matrix')
                imgData = obj.memmap_matrix_data(:,:,(frame-1)*obj.n_ch+channel);    
            else
                imgData = obj.memmap_matrix_data(:, (1:x_len/obj.n_ch) + (channel-1) * x_len/obj.n_ch, frame);
            end
            end

            %if strcmp(obj.type, 'binary') || strcmp(obj.type, 'matrix')
            %     imgData = imgData';
            %end
            if ~strcmp(obj.type, 'binary') %If object type is NOT binary, do not transpose(???)
                imgData = imgData;%';
            end

            if opt == 0
                set(data.im{channel}, 'CData', imgData);
            else
                imagesc(obj.ax{channel}, imgData);
                set(obj.ax{channel}, 'XTick', [], 'YTick', []);
                axis(obj.ax{channel}, 'image'); % Maintain aspect ratio
                data.im{channel}=obj.ax{channel}.Children;
                colormap('gray');
            end
            obj.image_object{channel}=obj.ax{channel}.Children;
            guidata(obj.figure,data);
        end
        
        function updateFrameFromFile(obj, channel, opt)
            data=guidata(obj.figure);
            obj.memmap_data = (obj.CurrFrame - 1) * obj.n_ch + channel;
            if ~iscell(obj.filename)
                datavals = bigread4(obj.filename, obj.memmap_data, obj.n_ch);
            else
                datavals = bigread4(obj.filename{obj.memmap_data}, 1, obj.n_ch);
            end
            imgData = datavals(:,:,1);

            if opt == 0
                set(data.im{channel}, 'CData', imgData);
            else
                imagesc(obj.ax{channel}, imgData);
                data.im{channel}=obj.ax{channel}.Children;
                set(obj.ax{channel}, 'XTick', [], 'YTick', []);
                axis(obj.ax{channel}, 'image'); % Maintain aspect ratio
            end
            guidata(obj.figure,data);
        end

function tv = setupTiffMapping(tv, filename, n_ch, info)
    % setupTiffMapping - Set up memory mapping for TIFF files
    %
    % Arguments:
    %   filename - Path to the TIFF file or folder containing TIFFs
    %   n_ch     - Number of channels
    %
    % Returns:
    %   tv - Struct containing memory-mapped TIFF data and related info
    
    % Check if filename is a string (single TIFF file)
    if ischar(filename)
        [folder, file, ext] = fileparts(filename);
        if isempty(ext)  % Case when we provide a folder, not a file
            % Read all TIFF files in the folder
            tv.memmap_data = tv.filename;
        else
            tv.filename = filename;
            if nargin<3 || isempty(info)
            info = readtifftags(filename);
            end
            offset_field = get_offset_field(info);
            
            % Determine if the TIFF frames have uniform offsets
            if length(info) > 3
                uneven_flag = info(2).(offset_field)(1) - info(1).(offset_field)(1) ~= info(3).(offset_field)(1) - info(2).(offset_field)(1);
            else
                uneven_flag = false;
            end
            
            if ~uneven_flag
                tv.map_type = 'mem';
                tv.memmap = memory_map_tiff(filename, [], n_ch, true);
                if isfield(info, 'GapBetweenImages') && info(1).GapBetweenImages == 0

                    % tv.memmap = memory_map_tiff(filename,[],n_ch,true);
                    tv.memmap_matrix = memory_map_tiff(filename,'matrix',n_ch,true);
                    tv.memmap_matrix_data = tv.memmap_matrix.Data.allchans;           
                end
                
                tv.width = info(1).ImageWidth;
                tv.height = info(1).ImageHeight;
            else
                tv.map_type = 'file';
            end
        end
    end
end

function tv = setupBinaryMapping(tv, filename, n_ch, framesize, form)
    % setupBinaryMapping - Set up memory mapping for binary files
    %
    % Arguments:
    %   filename  - Path to the binary file
    %   n_ch      - Number of channels
    %   framesize - Size of each image frame as [height, width]
    %   form      - Data format in the binary file (e.g., 'double', 'uint16')
    %
    % Returns:
    %   tv - Struct containing memory-mapped binary data and related info

    tv.filename = filename;
    tv.map_type = 'mem';
    
    % Prepare memory mapping format for each channel
    % format_string = cell(n_ch, 3);
    % for ch_rep = 1:n_ch
    %     format_string(ch_rep, :) = {form, framesize, ['channel', num2str(ch_rep)]};
    % end

    % Set up memory mapping

                if isempty(tv.n_ch)
                    userinputs = inputdlg({'Number of channels:','Height and Width (as matrix):','Data Format:'});
                    tv.n_ch = str2double(userinputs{1});
                    framesize = str2num(userinputs{2});
                    tv.height=framesize(1);
                    tv.width=framesize(2);
                    form=userinputs{3};
                end
    tv.memmap_matrix = memmapfile(filename, 'Format', {form, [framesize], 'allchans'}, 'Writable', false);
    n=length(tv.memmap_matrix.Data);
    tv.memmap_matrix = memmapfile(filename, 'Format', {form, [framesize n], 'allchans'}, 'Writable', false);
    % Store mapped data
    tv.memmap_matrix_data = tv.memmap_matrix.Data.allchans;    
    % Number of frames is determined by the binary file size
    tv.numFrames = n / tv.n_ch;
end

    end
end

function makeplot(hObject,event,tv)
data=guidata(tv.figure);
tv.CurrFrame=round(get(hObject,'Value'));
set(data.h.edit,'String',num2str(tv.CurrFrame));
guidata(tv.figure,data);
%displayFrame(tv);
tv.displayFrame(); % Changed to use object method
end

function makeplot2(hObject,event,tv)
data=guidata(tv.figure);
curval=get(hObject,'String');
try
    obj.CurrFrame=max(min(round(str2double(get(hObject,'String'))),tv.numFrames),1);
    hObject.String=num2str(obj.CurrFrame);
catch
    hObject.String=curval;
    return;
end
set(data.h.slide,'Value',obj.CurrFrame);
guidata(tv.figure,data);
%displayFrame(tv);
tv.displayFrame(); % Changed to use object method
end

function fps_but_down(hObject,event,tv)
data=guidata(tv.figure);
answer=inputdlg('Input FPS for Video Playback');
tv.fps=str2double(answer{1});
guidata(tv.figure,data);

end

function allimages = mm_proj(~, ~, tv, type)
    color_name = {'Red', 'Green', 'Blue'};
    f_out = figure;
    max_time = tv.max_time;
    
    % Set cursor to 'watch' during processing
    set([tv.figure, f_out], 'pointer', 'watch');
    drawnow;
    
    % Determine the number of subplots
    n_subplots = tv.n_ch + (tv.n_ch > 1); % Add extra plot for merged channels if more than one channel
    
    sub_handle_popup = zeros(1, n_subplots); % Initialize subplots
    P = cell(1, tv.n_ch); % Cell to store each channel's projection

    switch tv.map_type
        case 'mem'
            for ch = 1:tv.n_ch
                tic;
                sub_handle_popup(ch) = subplot(1, n_subplots, ch);
                
                % Compute projection based on type ('mean' or 'max')
                P{ch} = compute_projection(tv, type, ch, max_time, f_out, sub_handle_popup, n_subplots);
                
                % Display image for each channel
                if tv.n_ch==1
                    color_name{ch}=[];
                end
                display_image(P{ch}, f_out, sub_handle_popup(ch), tv.type, color_name{ch});
            end
            
            % Handle merged channel view if more than one channel
            if tv.n_ch > 1
                sub_handle_popup(end) = subplot(1, n_subplots, n_subplots);
                allimages = merge_channels(P, tv.type);
                imagesc(allimages); axis off;
                title('Merge');
            else
                allimages = P{1}';
            end
            
            disp('Projection Done.');
            linkaxes(sub_handle_popup); % Link the subplots for panning and zooming together
            
        case 'file'
            warning('Could not memory map, so projection would take too long.');
    end
    
    % Reset cursor to default
    set([tv.figure, f_out], 'pointer', 'arrow');
end

% Function to compute projection for each channel
function projection = compute_projection(tv, type, ch, max_time, f_out, sub_handle_popup, n_subplots)
    projection = [];
    
    if isempty(tv.memmap_matrix_data)
        projection = process_memmap(tv, ch, type, max_time, f_out, sub_handle_popup, n_subplots, 'memmap_data');
    else
        projection = process_memmap(tv, ch, type, max_time, f_out, sub_handle_popup, n_subplots, 'memmap_matrix_data');
    end
end

% Function to process memory-mapped data for projections
function P = process_memmap(tv, ch, type, max_time, f_out, sub_handle_popup, n_subplots, data_type)
    % Combined function to process memory-mapped data or matrix data based on type
    %
    % Arguments:
    %   tv                - Struct containing memory-mapped data
    %   ch                - Current channel being processed
    %   type              - Projection type (e.g., 'mean', 'max')
    %   max_time          - Maximum time allowed for processing
    %   f_out             - Figure handle for displaying images
    %   sub_handle_popup  - Subplot handles for displaying images
    %   n_subplots        - Number of subplots to display images
    %   data_type         - 'memmap_data' or 'memmap_matrix_data'
    
    tic;
    
    if strcmp(data_type, 'memmap_data')
        % Process memory-mapped data
        P = double(tv.memmap_data(1).(['channel', num2str(ch)])) / tv.numFrames;
        
        for b = 2:length(tv.memmap_data)
            P = imadd(P, double(tv.memmap_data(b).(['channel', num2str(ch)])) / tv.numFrames);
            
            if mod(b, 1000) == 0
                figure(f_out); subplot(1, n_subplots, ch);
                imagesc(P'); axis off; colormap('gray'); drawnow;
            end
            
            if toc > max_time / 2
                disp(['Max time reached for channel ', num2str(ch), '.']);
                P = P * (tv.numFrames / b);
                disp(['Frames averaged: ', num2str(b)]);
                break;
            end
        end
        
    elseif strcmp(data_type, 'memmap_matrix_data')
        % Process memory-mapped matrix data
        [y_len, x_len] = size(tv.memmap_matrix_data, 1:2);
        subdiv = 5; % Process in chunks to avoid memory overload
        n_subdiv = ceil(tv.numFrames / subdiv);
        P = zeros(x_len, y_len / tv.n_ch, 'double');
        
        for rep = 1:n_subdiv
            frames = 1 + (rep - 1) * subdiv : min(subdiv * rep, tv.numFrames);
            temp_frames = sum(tv.memmap_matrix_data(:, (1:x_len/tv.n_ch) + (ch-1) * x_len/tv.n_ch, frames), 3) / tv.numFrames;
            P = imadd(P, permute(temp_frames,[2 1 3]));

            if toc > max_time / 2
                disp(['Max time reached for channel ', num2str(ch), '.']);
                P = P / (rep * subdiv);
                disp(['Frames averaged: ', num2str(rep * subdiv)]);
                break;
            end
            
            figure(f_out); subplot(1, n_subplots, ch);
            imagesc(P'); axis off; colormap('gray'); drawnow;
        end
    else
        error('Unknown data_type specified. Use ''memmap_data'' or ''memmap_matrix_data''.');
    end
end

% Function to display an image
function display_image(image_data, f_out, subplot_handle, img_type, title_name)
    figure(f_out); subplot(subplot_handle);
    if strcmp(img_type, 'binary')
        imagesc(image_data); axis off; colormap('gray'); drawnow;
    else
        imagesc(image_data'); axis off; colormap('gray'); drawnow;
    end
    title(title_name);
end

% Function to merge channels for visualization
function merged_image = merge_channels(P, img_type)
    merged_image = permute(cat(3, P{:}), [2, 1, 3]);
    merged_image = merged_image ./ max(merged_image, [], [1, 2]);
    
    if size(merged_image, 3) == 2
        % Add a third zero-filled channel for RGB display if only 2 channels
        merged_image = cat(3, merged_image, zeros(size(merged_image, 1, 2), class(merged_image)));
    end
end

function offset_field=get_offset_field(info)
if isfield(info,'StripOffsets')
    offset_field='StripOffsets';
elseif isfield(info,'TileOffsets')
    offset_field='TileOffsets';
else
    error('Neither strip nor tile format.')
end
end

function play_vid(hObject,event,tv)
data=guidata(tv.figure);
frame_rate=1./hObject.InstantPeriod;
data.increment=max(floor(tv.fps/frame_rate),1);
tv.CurrFrame=mod(tv.CurrFrame+data.increment,tv.numFrames);
if tv.CurrFrame==0
    tv.CurrFrame=tv.numFrames;
end
data.h.slide.Value=tv.CurrFrame;
data.h.edit.String=num2str(tv.CurrFrame);
displayFrame(tv);
guidata(tv.figure,data);
end

function play_but_down(hObject,event,tv)
data=guidata(tv.figure);
if tv.CurrFrame==tv.numFrames
    tv.CurrFrame=1;
end
guidata(tv.figure,data);
set(data.timer,'Period',max(round(1/tv.fps,3),.001),'TimerFcn',{@play_vid,tv});
start(data.timer);
guidata(tv.figure,data);
set(hObject,'callback',@(x,evt) stop_but_down(x,evt,tv));
set(hObject,'String','=');
end

function stop_but_down(hObject,event,tv)
data=guidata(tv.figure);
stop(data.timer);
set(hObject,'callback',@(x,evt) play_but_down(x,evt,tv));
set(hObject,'String','>');
guidata(tv.figure,data);
end