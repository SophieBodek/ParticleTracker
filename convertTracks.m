function vtracks_convert = convertTracks(track_data, fps, ppx, ppy)
% Usage: vtracks_convert = convertTracks(track_data, fps, ppx, ppy)
% Converts units for vtracks struct outputted by PredictiveTracker function.
% Must input struct with fields: len, X, Y, U, V, T (will also accept
% struct with vtracks field containing len, X, Y, U, V, T). Also requires
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

% Check that track_data has correct structure (fields: len, X, Y, U, V, T)
if isfield(track_data, 'vtracks') 
    % if track_data has a vtracks field
    vtracks = track_data.vtracks;
    fprintf('Input format: struct with vtracks field\n');
elseif isfield(track_data, 'X') && isfield(track_data, 'Y')
    % if track_data is directly the vtracks struct
    vtracks = track_data;
    fprintf('Input format: direct vtracks struct\n');
end

% Check that vtracks has correct fields
required_fields = {'len', 'X', 'Y', 'U', 'V', 'T'};
struct_fields = fieldnames(vtracks);
if ~all(ismember(required_fields, struct_fields))
    missing_fields = required_fields(~ismember(required_fields, struct_fields));
    error('Invalid vtracks struct. Missing required fields: %s', strjoin(missing_fields, ', '));
end
    
vtracks_convert = vtracks; % Initialize output struct

ntracks = numel(vtracks);
fprintf('%d tracks to convert\n', ntracks);

for tr = 1:ntracks
    vtracks_convert(tr).X = vtracks(tr).X ./ ppx; %convert from pixels to units x direction
    vtracks_convert(tr).Y = vtracks(tr).Y ./ ppy; %convert from pixels to units y direction
    vtracks_convert(tr).U = vtracks(tr).U .* (fps/ppx); %convert velocity to units per second
    vtracks_convert(tr).V = vtracks(tr).V .* (fps/ppx); 
    vtracks_convert(tr).T = vtracks(tr).T ./ fps; %convert frames to seconds

    if mod(tr, 10) == 0 || tr == ntracks
        fprintf('Converted %d/%d tracks\n', tr, ntracks);
    end
end

fprintf('All %d tracks converted!\n', ntracks);   
end