function [vtracks_filt, ntracks_filt] = FilterTracks(vtracks, lenMin, xDistMin, yDistMin, lenComp, xComp, yComp, xbounds, ybounds, tbounds, split_behavior, noisy)
% FilterTracks - Filters particle tracks based on length and displacement criteria;
%                can select for tracks 'greater' or 'less' than a specified
%                track length (lenMin), greater/less than a specified X or
%                Y displacement (xDistMin, yDistMin), as well as tracks within 
%                specified spatial or time boundaries (xbounds, ybounds, tbounds). 
%                Split_behavior indicates treatment of tracks that start within, 
%                but finish outside of given boundaries ('include', 'exclude', or
%                'trim'). Can generate plot of filtered tracks (noisy = 1).
%
% Inputs:
%   vtracks   - Input track data structure
%   lenMin    - Minimum track length threshold (default: 20); track length
%               refers to number of data points in the track
%   xDistMin  - X displacement threshold (default: 0)
%   yDistMin  - Y displacement threshold (default: 0)
%   lenComp   - Length comparison: 'greater' (default) or 'less'
%   xComp     - X displacement comparison: 'greater' (default) or 'less'
%   yComp     - Y displacement comparison: 'greater' (default) or 'less'
%   xbounds   - X coordinate boundaries [min, max] (empty to ignore)
%   ybounds   - Y coordinate boundaries [min, max] (empty to ignore)
%   tbounds   - T spatial coordinate limits [start_time, end_time] 
%               (empty to ignore)
%   split_behavior - Boundary handling: 'exclude' (default), 'include', or 'trim'
%                    'exclude': reject tracks with any points outside bounds
%                    'include': accept tracks that start within bounds
%                    'trim': truncate tracks to only points within bounds
%   noisy     - Boolean for plotting results

% Default values
lenMin_default = 20; % typically number of frames
xDistMin_default = 0; % distance units
yDistMin_default = 0; % distance units
noisy_default = 0; % no plotting

% Set default values - generates printout if defaults are used
if nargin < 1
    error(['Usage: [vtracks_filt, ntrack_filt] = ' mfilename ...
        '(vtracks, lenMin, xDistMin, yDistMin, lenComp, xComp, yComp, xbounds, ' ...
        'ybounds, tbounds, split_behavior, noisy)']);
end
if nargin < 2 || isempty(lenMin)
    lenMin = lenMin_default; 
    fprintf('No track length threshold specified; using length = %d\n', lenMin_default);
end
if nargin < 3 || isempty(xDistMin)
    xDistMin = xDistMin_default;
    fprintf('No X displacement threshold specified; using %d [units]\n', xDistMin_default);
end
if nargin < 4 || isempty(yDistMin)
    yDistMin = yDistMin_default; 
    fprintf('No Y displacement threshold specified; using %d [units]\n', yDistMin_default);
end
if nargin < 5 || isempty(lenComp), lenComp = 'greater'; end %default is greater than given length
if nargin < 6 || isempty(xComp), xComp = 'greater'; end %default is greater than given X displacement
if nargin < 7 || isempty(yComp), yComp = 'greater'; end %default is greater than given Y displacement
% 8, 9, 10 are xbounds, ybounds, tbounds
if nargin < 11 || isempty(split_behavior), split_behavior = 'exclude'; end % default is exclude out of bounds tracks
if ~exist('noisy','var') || isempty(noisy)
        noisy = noisy_default;
end

% Check that track input has correct structure (fields: len, X, Y, U, V, T)
if isfield(vtracks, 'vtracks') 
    vtracks = vtracks.vtracks; %track data may be in a vtracks field
elseif isfield(vtracks, 'len') && isfield(vtracks, 'X')
    vtracks = vtracks; %otherwise vtracks should have the fields we are looking for
end
% Check that track input has ALL the correct fields (at least len,X,Y,T)
required_fields = {'len','X', 'Y','T'};
struct_fields = fieldnames(vtracks);
if ~all(ismember(required_fields, struct_fields))
    missing_fields = required_fields(~ismember(required_fields, struct_fields));
    error('Invalid vtracks struct. Missing required fields: %s', strjoin(missing_fields, ', '));
end
all_fields = fieldnames(vtracks); % preserves U, V, Ax, Ay, Theta (and more) if included

% Validate comparison inputs
valid_comparisons = {'greater', 'less'};
if ~ismember(lenComp, valid_comparisons)
    error('FilterTracks:InvalidLenComp', ...
          'lenComp must be ''greater'' or ''less''. Got: %s', lenComp);
end
if ~ismember(xComp, valid_comparisons)
    error('FilterTracks:InvalidXComp', ...
          'xComp must be ''greater'' or ''less''. Got: %s', xComp);
end
if ~ismember(yComp, valid_comparisons)
    error('FilterTracks:InvalidYComp', ...
          'yComp must be ''greater'' or ''less''. Got: %s', yComp);
end

% Check split_behavior has appropriate keyword
valid_behaviors = {'exclude', 'include', 'trim'};
if ~ismember(split_behavior, valid_behaviors)
    error('FilterTracks:InvalidSplitBehavior', ...
          'split_behavior must be ''exclude'', ''include'', or ''trim''. Got: %s', split_behavior);
end

% FILTERING TRACKS --------------------------------------------------------
% -------------------------------------------------------------------------
disp('Filtering tracks...');

vtracks_filt = [];
trackLens = extractfield(vtracks,'len');
count = 0;

% Filter by track length --------------------------------------------------
if strcmp(lenComp, 'greater') %greater than or equal to specified length
    validLenInd = find(trackLens >= lenMin); %done to match with longtracks
    fprintf('- finding tracks greater than length = %d\n', lenMin)
    if isempty(validLenInd) %warning of no tracks found
        warning('FilterTracks:NoLengthMatch', ...
                'No tracks found with length greater than %d. Maximum track length is %d.', ...
                lenMin, max(trackLens));
        %return;
    end
else
    validLenInd = find(trackLens < lenMin);
    fprintf('- finding tracks less than length = %d\n', lenMin)
    if isempty(validLenInd)
        warning('FilterTracks:NoLengthMatch', ...
                'No tracks found with length less than %d. Minimum track length is %d.', ...
                lenMin, min(trackLens));
        %return;
    end
end

% Filter tracks by spatial and time boundaries ----------------------------
initialCount = length(validLenInd);
xBoundRejects = 0;
yBoundRejects = 0;
tBoundRejects = 0;
xDispRejects = 0;
yDispRejects = 0;
trimmedTracks = 0;

for i = 1:length(validLenInd)
    trackInd = validLenInd(i);
    x = vtracks(trackInd).X;
    y = vtracks(trackInd).Y;
    t = vtracks(trackInd).T;
    
    valid_indices = true(size(x)); % all points initially valid
    
    switch split_behavior
        case 'exclude' %exclude tracks if 1+ point(s) are outside bounds
            if ~isempty(xbounds) %X boundaries
                if sum(x > max(xbounds)) > 1 || sum(x < min(xbounds)) > 1
                    xBoundRejects = xBoundRejects + 1;
                    continue
                end
            end
            if ~isempty(ybounds) %Y boundaries
                if sum(y > max(ybounds)) > 1 || sum(y < min(ybounds)) > 1
                    yBoundRejects = yBoundRejects + 1;
                    continue
                end
            end
            if ~isempty(tbounds) %T boundaries
                if sum(t > max(tbounds)) > 1 || sum(t < min(tbounds)) > 1
                    tBoundRejects = tBoundRejects + 1;
                    continue
                end
            end
            
        case 'include' %include tracks if starts within bounds
            if ~isempty(xbounds) %X boundaries
                if x(1) < min(xbounds) || x(1) > max(xbounds)
                    xBoundRejects = xBoundRejects + 1;
                    continue
                end
            end
            if ~isempty(ybounds) %Y boundaries
                if y(1) < min(ybounds) || y(1) > max(ybounds)
                    yBoundRejects = yBoundRejects + 1;
                    continue
                end
            end
            if ~isempty(tbounds) %T boundaries
                if t(1) < min(tbounds) || t(1) > max(tbounds)
                    tBoundRejects = tBoundRejects + 1;
                    continue
                end
            end
            
        case 'trim' %trim tracks that cross boundaries
            if ~isempty(xbounds)
                valid_indices = valid_indices & (x >= min(xbounds)) & (x <= max(xbounds));
            end
            if ~isempty(ybounds)
                valid_indices = valid_indices & (y >= min(ybounds)) & (y <= max(ybounds));
            end
            if ~isempty(tbounds)
                valid_indices = valid_indices & (t >= min(tbounds)) & (t <= max(tbounds));
            end
            
            if sum(valid_indices) == 0
                xBoundRejects = xBoundRejects + 1; % Count as boundary reject
                continue
            end
            
            x = x(valid_indices); %trimming
            y = y(valid_indices);
            t = t(valid_indices);
            
            if sum(valid_indices) < length(valid_indices)
                trimmedTracks = trimmedTracks + 1;
            end
    end
    
    % Calculate displacements - based on trimmed tracks
    x1 = x(1);
    x2 = x(end);
    y1 = y(1);
    y2 = y(end);
    
    xDist = abs(x2 - x1);  % Absolute X displacement
    yDist = abs(y2 - y1);  % Absolute Y displacement
    
    % Filter by X, Y displacements ------------------
    if strcmp(xComp, 'greater')
        if xDist < xDistMin
            xDispRejects = xDispRejects + 1;
            continue
        end
    else
        if xDist > xDistMin
            xDispRejects = xDispRejects + 1;
            continue
        end
    end
    if strcmp(yComp, 'greater')
        if yDist < yDistMin
            yDispRejects = yDispRejects + 1;
            continue
        end
    else
        if yDist > yDistMin
            yDispRejects = yDispRejects + 1;
            continue
        end
    end
    
    %Total track counter
    count = count + 1;
    
    if strcmp(split_behavior, 'trim') %use trimmed data
        for fn = 1:length(all_fields)
            field = all_fields{fn};
        
            if strcmp(field, 'len')
                % Special handling for length
                vtracks_filt(count).len = sum(valid_indices);
            else
                % Trim all other fields using valid_indices
                vtracks_filt(count).(field) = vtracks(trackInd).(field)(valid_indices);
            end
        end
    else % include/exclude - all fields from original track
        for fn = 1:length(all_fields)
            field = all_fields{fn};
            vtracks_filt(count).(field) = vtracks(trackInd).(field);
        end
    end
end

% Display outcomes and issue warnings
ntracks_filt = length(vtracks_filt);

if ~isempty(xbounds), fprintf('- finding tracks within X boundaries [%.1f, %.1f]\n', min(xbounds), max(xbounds)); end
if ~isempty(xbounds) && xBoundRejects == initialCount
    warning('FilterTracks:NoXBoundaryMatch', ...
        'No tracks remain within X boundaries [%.1f, %.1f]. All %d tracks were rejected.', ...
        min(xbounds), max(xbounds), xBoundRejects);
elseif xBoundRejects > 0
    warning('FilterTracks:SomeXBoundaryRejects', ...
        '%d out of %d tracks were rejected for exceeding X boundaries [%.1f, %.1f].', ...
        xBoundRejects, initialCount, min(xbounds), max(xbounds));
end

if ~isempty(ybounds), fprintf('- finding tracks within Y boundaries [%.1f, %.1f]\n', min(ybounds), max(ybounds)); end
if ~isempty(ybounds) && yBoundRejects == initialCount - xBoundRejects
    warning('FilterTracks:NoYBoundaryMatch', ...
        'No tracks remain within Y boundaries [%.1f, %.1f]. All remaining tracks were rejected.', ...
        min(ybounds), max(ybounds));
elseif yBoundRejects > 0
    warning('FilterTracks:SomeYBoundaryRejects', ...
        '%d tracks were rejected for exceeding Y boundaries [%.1f, %.1f].', ...
        yBoundRejects, min(ybounds), max(ybounds));
end

if ~isempty(tbounds), fprintf('- finding tracks within T boundaries [%.1f, %.1f]\n', min(tbounds), max(tbounds)); end
if ~isempty(tbounds) && tBoundRejects == initialCount - xBoundRejects - yBoundRejects
    warning('FilterTracks:NoTBoundaryMatch', ...
        'No tracks remain within T boundaries [%.1f, %.1f]. All remaining tracks were rejected.', ...
        min(tbounds), max(tbounds));
elseif tBoundRejects > 0
    warning('FilterTracks:SomeTBoundaryRejects', ...
        '%d tracks were rejected for exceeding T boundaries [%.1f, %.1f].', ...
        tBoundRejects, min(tbounds), max(tbounds));
end

remaining = initialCount - xBoundRejects - yBoundRejects - tBoundRejects;
fprintf('- finding tracks with X displacement %s than %d [units]\n', xComp, xDistMin)
if xDispRejects == remaining
    warning('FilterTracks:NoXDisplacementMatch', ...
        'No tracks have X displacement %s than %.1f. All remaining tracks were rejected.', ...
        xComp, xDistMin);
elseif xDispRejects > 0
    warning('FilterTracks:SomeXDisplacementRejects', ...
        '%d tracks were rejected for X displacement not %s than %.1f.', ...
        xDispRejects, xComp, xDistMin);
end

remaining = remaining - xDispRejects;
fprintf('- finding tracks with Y displacement %s than %d [units]\n', yComp, yDistMin)
if yDispRejects == remaining
    warning('FilterTracks:NoYDisplacementMatch', ...
        'No tracks have Y displacement %s than %.1f. All remaining tracks were rejected.', ...
        yComp, yDistMin);
elseif yDispRejects > 0
    warning('FilterTracks:SomeYDisplacementRejects', ...
        '%d tracks were rejected for Y displacement not %s than %.1f.', ...
        yDispRejects, yComp, yDistMin);
end

% Final summary warning if no tracks pass all filters
if count == 0
    warning('FilterTracks:NoTracksRemain', ...
        'No tracks passed all filtering criteria. Started with %d tracks, all were filtered out.', ...
        length(vtracks));
else
    fprintf('%d tracks fall within given parameters. \n', ntracks_filt);
    if strcmp(split_behavior, 'trim') && trimmedTracks > 0
        fprintf('%d tracks were trimmed to fit within boundaries.\n', trimmedTracks);
    end
end

% Plotting if requested--------------------------------------------------
if noisy
    figure;
    hold on;
    for i = 1:ntracks_filt
        plot(vtracks_filt(i).X, vtracks_filt(i).Y,'LineWidth',1);
    end

    axis tight
    %set(gca, 'YDir','reverse') % reverse y-axis for better comparison with images (which have 0,0 in upper left)
    set(gca,'dataaspectratio',[1 1 1]) %preserve aspect ratio of images
    
    % Update title to reflect split_behavior
    fig_title = {['Track lengths ' lenComp ' than ' int2str(lenMin) '; X displacement '],
                 [xComp ' than ' int2str(xDistMin) ' [units]; Y displacement ' yComp ...
                  ' than ' int2str(yDistMin) ' [units].']};
    title(fig_title,'interpreter','none')
end
end