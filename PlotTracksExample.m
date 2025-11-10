% PlotTracksExample.m provides code to document support functions for
% PredictiveTracker.m, as well as demonstrate built-in visualization
% tools for examining the particle track data structure output of
% DifferentiateTracks.m (or PredictiveTracker.m) functions.
%
% This code requires/discusses the following functions:
% -=--- XXX -=------------------------------------------------------------
%   -   BackgroundImage : generates background image for ParticleFinder
%   -   ParticleFinder : identifies particle locations; also plays a video 
%                        of particle locations.
%   -   PredictiveTracker : links particle coordinates into tracks
%                           (ParticleFinder and DifferentiateTracks are
%                           built in; see TrackExample.m); can play a video
%                           of particle trajectories.
%   -   DifferentiateTracks : calculates velocities and accelerations from
%                             PredictiveTracker trajectories; can also
%                             output smoothed trajectories. 
%   -   ConvertTracks : converts tracks output from pixels and frames into
%                       specified spatial and time units.
%   -   FilterTracks : filters tracks based on duration, displacement, or
%                      location; also plots filtered tracks.
%   -   longtracks.m : plots tracks over a given length (i.e., number of frames)
%                      with two plot types 
% -=--- XXX -=------------------------------------------------------------
%   -   Velocities : outputs [X, Y, T, U, V, Tr] from vtracks
%                    also plots particle velocity vectors
%   -   Accelerations : outputs [Ax, Ay, X, Y, T, Tr] from vtracks; 
%                       also plots particle accelerations
%   -   Curvature : outputs [Curv, X, Y, T, Tr] from vtracks
%   -   Vorticities : outputs [Vort, X, Y, T, Tr] from vtracks
%                     also plots particle vorticities
%                     (requires Partial Differential Equation toolbox)
%   -   Divergences : outputs [Div, X, Y, T, Tr] from vtracks
%                     also plots particle divergences
%                     (requires Partial Differential Equation toolbox)
% -=--- XXX -=------------------------------------------------------------
%   -   TiffTrackViewer : viewer for playing images as a movie with tracks 
%                         overlaid; has functionality to scroll back/forward
%                         through the movie and to click on tracks for 
%                         individual trajectory plotting.
%   -   bigread4 : quickly loads large .TIFF stacks (required for TiffTrackViewer)


% -=- Make a background image -=------------------------------------------
%     Makes an ensenble average of all frames so we only track moving objects
im_file = 'TrackExample/0*.tif'; % directory with images
bg_file = 'TrackExample/background.tif'; % name of background file

%bg = BackgroundImage(im_file, bg_file);


% -=- Run FindParticles to identify particle coodinates -=----------------
%     Identifies particles in a series of images (im_file) based on a brightness
%     threshold (threshold) above the background (bg_file); you can also specify
%     a range of frames to for finding particles (framerange), approximate area 
%     of particles (arealim), whether to search for particles that are 
%     brighter/darker than the background (invert), and whether to play a video of 
%     particle locations (noisy).
threshold = 12; 
framerange = [1 inf];
arealim = 1; % single pixel area
invert = 0; % brighter than background
noisy = 0; %plotting/video parameter : no plot = 0; plot or play video= 1

[x,y,t] = ParticleFinder(im_file, threshold, framerange, [], bg_file, arealim, invert, noisy);
% TO DO: fix output to file; check function still works for .avi filetype;
%        fix saving output video frames (only saves lower corner--see PredTrackerMovie.m)


% -=- Use PredictiveTracker to link particle coordinates into tracks -=----
%     PredictiveTracker can be used in two ways: (1) you can input a series
%     of images with specified particle-identifying parameters (same as for
%     FindParticles above) and particle-tracking parameters (max_disp
%     indicates maximum search radius in predicted position); (2) you can
%     use a separate particle-identifying algorithm and feed the particle
%     coordinates into PredictiveTracker in an array of the form [x y t],
%     then use the same particle-tracking parameters. Additionally,
%     PredictiveTracker can output either raw tracks (difftracks = 0), or 
%     perform built-in smoothing and differentiation to calculate particle 
%     velocities and accelerations depending on use-case. 
max_disp = 4;
minarea = arealim;
noisy = 0; 
difftracks = 0; % output raw tracks only

%    1.) find particles directly in PredictiveTracker
%[tracks,ntracks,lmean,lrms] = PredictiveTracker(im_file, max_disp, threshold, ...
%                                bg_file, minarea, invert, noisy, difftracks);

%    2.) input particle coordinates for tracking
[tracks,ntracks,lmean,lrms] = PredictiveTracker([x y t], max_disp, [], [], [], [], noisy, difftracks);
% TO DO: check function still works for .avi filetype; fix saving output video frames 
%        (only saves lower corner--see PredTrackerMovie.m); 
%        add movie/saving frames when not reading direct from images


% -=- Use DifferentiateTracks to calculate velocities and accelerations -=-
%     DifferentiateTracks is built into PredictiveTracker, but can be used
%     independently to perform smoothing/filtering of trajectories
%     themselves, in addition to smoothing and differentiation for
%     velocities and accelerations.
filterwidth = 2; % kernel width
fitwidth = 3; % fit window size
smoothed = 1; % output smoothed trajectories

vtracks = DifferentiateTracks(tracks, filterwidth, fitwidth, smoothed);


% -=- Use FilterTracks to select a specific subset of tracks -=-----------
%     We can filter trajectories based on duration (lenMin), horizontal and
%     vertical displacement (xDistMin, yDistMin), physical location in the
%     images (xbounds, ybounds) and a temporal subset of images (tbounds); 
%     boundary handling can also be specified (split_behavior). Further, we 
%     can indicate whether we want track durations or displacements to be 
%     "greater" or "less" than the specified limits (lenComp, xComp, yComp);
%     i.e, we can specify we want tracks that last longer than 100 frames,
%     but whose displacement does not exceed 5 distance units.
%     We can also plot the filtered tracks (noisy)
lenMin = 10;
xDistMin = 0; 
yDistMin = 0; % default is "greater" than lenMin, xDistMin and yDistMin
xbounds = [0 180]; % select tracks in bottom left corner
ybounds = [0 110];
tbounds = [1 inf]; % all frames
split_behavior = "trim"; % trim tracks outside the specified boundaries
noisy = 1; 

ftracks = FilterTracks(vtracks,lenMin, xDistMin, yDistMin, [], [], [], ...
                       xbounds, ybounds, tbounds, split_behavior, noisy);


% -=- Use longtracks to view tracks over minlength -=----------------------
%     Similar to the plotting output of FilterTracks, longtracks plots all
%     tracks longer than minlength. Additionally, when the plotting
%     parameter noisy = 2, the trajectories are shifted to begin at the
%     origin. This can help with visualizing diffusion or other transport
%     patterns.
minlength = 20; %minimum track length to plot
noisy = 2; % shifts tracks to origin

[nlong,ntracks,lmean,lrms,lmax]=longtracks(vtracks, minlength, noisy);


% -=- Convert units using convertTracks -=--------------------------------
%     ParticleFinder and PredictiveTracker output particle locations in
%     pixels and frames. ConvertTracks is a convenient way to convert your
%     tracks/vtracks struct into physical units. ConvertTracks can handle
%     inputs with the following quantities: time (T), position (X and Y), 
%     velocity (U and V), acceleration (Ax and Ay), curvature (K). 
fps = 30; %frames per second
ppx = 87; %pixels per unit (for example, 50 pixels per cm)
ppy = 87;

ctracks = convertTracks(vtracks, fps, ppx, ppy);


% Use built-in plotting functionality of Velocities, Accelerations, Vorticities,  
% and Divergences to visualize particle trajectories -----------------------------
framerange=[1 inf]; % all frames
noisy = 1; % show plots

% -=- Velocities -=-------------------------------------------------------
%     The Velocities function sorts the data in vtracks (x, y, u, v, t, tr)
%     by time and outputs each in a vector of the same name; can also generate
%     a quiver plot with particle velocities.

[u,v,x,y,t,tr]=Velocities(vtracks,framerange,noisy);


% -=- Accelerations -=----------------------------------------------------
%     The Accelerations function sorts the data in vtracks (x, y, ax, ay,
%     t, tr) by time and outputs each in a vector of the same name.
%     Also plots particle accelerations if requested.

[ax,ay,~,~,~,~]=Accelerations(vtracks,framerange,noisy);


% -=- Curvature -=-------------------------------------------------------
%     The Curvature function calculates curvature and outputs data from
%     vtracks into a data struct with additional field "K"; additionally 
%     outputs data into vectors that are sorted by time: k, x, y, t, tr. 
%     Also plots curvatures.

[ktracks,k,~,~,~,~]=Curvature(vtracks,framerange,noisy);


% -=- Vorticities -=------------------------------------------------------
%     The Vorticities function calculates vorticity and outputs data from
%     vtracks in vectors that are sorted by time: vort, x, y, t, tr.
%     Also plots particle vorticities. 

[vort,~,~,~,~]=Vorticities(vtracks,framerange,noisy); 


% -=- Divergences -=------------------------------------------------------
%     The Divergences function calculates divergence and outputs data from
%     vtracks in vectors that are sorted by time: div, x, y, t, tr. Can
%     also plot particle divergences.

[div,~,~,~,~]=Divergences(vtracks,framerange,noisy); 


% -=- TiffTrackViewer -=--------------------------------------------------
%     The TiffTrackViewer is a Matlab class that plays Tiff images or Tiff
%     stacks as a movie with overlaid tracks. The viewer can scroll
%     backward/forward through the movie and has functionality to select
%     individual tracks for plotting. 

% Add field for track index (helps with referencing track number in original vtracks struct)
ids = num2cell(1:length(vtracks));
[vtracks.ID] = deal(ids{:});

% Are the images saved as a .TIFF stack? or as a series of .TIFF files?
isStack = 1; %0 for invidiual images in a directory; 1 for .TIFF stack; 2 for both

if isStack == 0 || isStack == 2     %--------- Individual .TIFF images
    im_path = 'TrackExample/0*.tif';
    im_list = dir(im_path);
    names = extractfield(im_list,'name');
end

if isStack ==1 || isStack ==2     %--------- .TIFF stack
    stack_name = 'TrackExampleStack.tiff';

    % Parameters for reading .TIFF stacks
    tiff_frame_start = 1; %starting frame of tiff stack
    tiff_num_frames = []; %number of tiff frames to read after starting index

    TS = bigread4(stack_name, tiff_frame_start, tiff_num_frames); %load .TIFF stack
end

% Displays viewer with intuitive controls for playing images as a movie and
% overlaying tracks, as well as plotting trajectory metrics. 
n_ch = 1; %number of channels (1 if greyscale; 3 if color)
minlen = 20;  %minimum track length; if too low, viewer gets very laggy

% TiffTrackViewer : TIFF movie viewer and interactive track visualizer
if isStack == 0 %functionality more limited for directory with images
    tv = TiffTrackViewer('TrackExample', n_ch, tracks, minlen);
else %otherwise should display .TIFF stack
    tv = TiffTrackViewer(TS,n_ch,tracks,minlen); 
end

%% Visualization Method : Predictive Tracker Movie
% This function plays the same movie as the PredictiveTracker function
%noisy = 1; %noisy = 2 saves indivual frames

%disp('PredTrackerMovie : Plays video with active tracks');
%if isStack == 0 || isStack == 2 %directory with images (allows for saving movie frames)
%    PredTrackerMovie(im_path,vtracks);
%else % .TIFF stack
%    PredTrackerMovie(stack_path,vtracks);
%end