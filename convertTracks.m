function vtracks_convert = convertTracks(track_data, fps, ppx, ppy)
% Usage: vtracks_convert = convertTracks(track_data, fps, ppx, ppy)
% Converts units for vtracks struct outputted by PredictiveTracker function.
% Must input struct with fields: len, X, Y, and T (will also accept struct
% with vtracks fields for U, V, Ax, Ay, K and Theta). Also requires
% fps - frames per second for time conversion; and pixels per unit in the X
% and Y directions (ppx and ppy). 

fps_default = 30; %default frame rate is 30fps 

if nargin < 1
    error(['Usage: vtracks_convert = ' mfilename ...
        '(track_data, fps, ppx, ppy)']);
end
if ~exist('fps','var') || isempty(fps)
    fps = fps_default;
    fprintf('No fps indicated; using default %d fps', fps_default);
end
if ~exist('ppx','var') || isempty(ppx)
        error('Missing conversion: pixels per unit X direction');
end
if ~exist('ppy','var') || isempty(ppy)
        error('Missing conversion: pixels per unit Y direction');
end

% Check that track_data has correct structure (fields: len, X, Y, T)
if isfield(track_data, 'vtracks') 
    % if track_data has a vtracks field
    vtracks = track_data.vtracks;
elseif isfield(track_data, 'X') && isfield(track_data, 'Y')
    % if track_data is directly the vtracks struct
    vtracks = track_data;
end

% Check that vtracks has correct fields (needs X, Y, T at minimum)
required_fields = {'len', 'X', 'Y', 'T'};
struct_fields = fieldnames(vtracks);
if ~all(ismember(required_fields, struct_fields))
    missing_fields = required_fields(~ismember(required_fields, struct_fields));
    error('Invalid vtracks struct. Missing required fields: %s', strjoin(missing_fields, ', '));
end

% Check for optional fields
has_U = isfield(vtracks, 'U'); %velocities
has_V = isfield(vtracks, 'V');
has_Ax = isfield(vtracks, 'Ax'); %accelerations
has_Ay = isfield(vtracks, 'Ay');
has_K = isfield(vtracks, 'K'); %curvature
has_Theta = isfield(vtracks, 'Theta'); %angle/theta
    
vtracks_convert = vtracks; % Initialize output struct

ntracks = numel(vtracks);
fprintf('%d tracks to convert\n', ntracks);

for tr = 1:ntracks
    % convert required fields
    vtracks_convert(tr).X = vtracks(tr).X ./ ppx; %convert from pixels to units x direction
    vtracks_convert(tr).Y = vtracks(tr).Y ./ ppy; %convert from pixels to units y direction
    vtracks_convert(tr).T = vtracks(tr).T ./ fps; %convert frames to seconds

    if has_U % convert if velocity fields included
        vtracks_convert(tr).U = vtracks(tr).U .* (fps/ppx); % pixels per frame to units per second
    end
    if has_V
        vtracks_convert(tr).V = vtracks(tr).V .* (fps/ppy); 
    end
    
    if has_Ax % convert if accelerations fields included
        vtracks_convert(tr).Ax = vtracks(tr).Ax .* (fps^2/ppx); % px/frame^2 to units/s^2
    end
    if has_Ay
        vtracks_convert(tr).Ay = vtracks(tr).Ay .* (fps^2/ppy); 
    end
    
    if has_K % convert curvature if included
        if ppx ~= ppy % curvature conversion is wonky if horizontal and vertical scales differ
            warning('Horizontal and vertical scale differ; consider converting tracks before calculating curvature.');
        end
        vtracks_convert(tr).K = vtracks(tr).K .* sqrt(ppx*ppy); % 1/pixel to 1/units
    end
    
    if has_Theta % preserve Theta (no conversion)
        vtracks_convert(tr).Theta = vtracks(tr).Theta;
    end

    if mod(tr, 10) == 0 || tr == ntracks
        fprintf('Converted %d/%d tracks\n', tr, ntracks);
    end
end

fprintf('All %d tracks converted!\n', ntracks);   
end