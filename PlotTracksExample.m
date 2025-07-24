% PlotTracks.m, written July 2025 by Sophie Bodek, provides code that can
% visualize the vtracks output of the PredictiveTracker.m function using
% different methods, such as the longtracks.m function or
% PredTrackerMovie.m function. 
% 
% In order to run, this does require the following:
%   -   bigread4.m : quickly loads large .TIFF stacks
%   -   Velocities.m : outputs X, Y, T, U, V, Tr from vtracks
%                      also plots particle velocity vectors
%   -   TiffTrackViewer.m : viewer for playing images as a movie with tracks overlaid
%   -   longtracks.m : plots tracks over a given length
%   -   PredTrackerMovie.m : plays active tracks over images as a movie --> TO ADD MAYBE

close all
clear all

%% %--------- FOLDERS AND FILENAMES
% Are the images saved as a .TIFF stack? or as a series of .TIFF files?
isStack = 2; %0 for invidiual images in a directory; 1 for .TIFF stack; 2 for both

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

%--------- Load vtracks 
% (assumes you have already run PredictiveTracker elsewhere)
load('vtracks.mat');

% Add field for track index (helps with referencing track number in original vtracks struct)
ids = num2cell(1:length(vtracks));
[vtracks.ID] = deal(ids{:});

%--------- Calculate velocities
framerange=[1 inf];
noisy = 1; %plotting parameter : no plot = 0; plot velocity vectors = 1
% Note: Velocities.m also has a visualization tool (noisy = 1) that plots
% the velocity vector of each patricle
disp('Velocities : plots particle velocity vectors');
[u,v,x,y,t,tr]=Velocities(vtracks,framerange,noisy);

%fig_velocities = gcf; %save velocities figure if desired
%saveas(fig_velocities,[directory filesep 'velocities.png']); 

disp('Tracks locked & loaded!')

%% Visualization Method 1: TIFF TRACK VIEWER
% Displays viewer with intuitive controls for playing images as a movie and
% overlaying tracks, as well as plotting trajectory metrics
n_ch = 1; %number of channels (1 if greyscale; 3 if color)
minlen = 50;  %minimum track length

disp('TiffTrackViewer : TIFF movie viewer and interactive track visualizer');
if isStack == 0 %functionality more limited for directory with images
    tv = TiffTrackViewer('TrackExample', n_ch, vtracks, minlen);
else %otherwise should display .TIFF stack
    tv = TiffTrackViewer(TS,n_ch,vtracks,minlen); 
end

%% Visualization Method 2: Long Tracks
% This function plots all tracks over a given length
tracks = sortrows([tr x y], 1); %need to sort data due to how function is written
noisy = 1; %display plot
minlength = 50; %minimum track length to plot
readwhole = 0; %'True';

fprintf('longtracks : Plots all tracks longer than %d.\n', minlength);
[nlong,ntracks,lmean,lrms,lmax]=longtracks(tracks,minlength,noisy);

%fig_longtracks = gcf; %to save the figure
%saveas(fig_longtracks,[directory filesep 'longtracks.png']); 

%% Visualization Method 3: Predictive Tracker Movie
% This function plays the same movie as the PredictiveTracker function
%noisy = 1; %noisy = 2 saves indivual frames

%disp('PredTrackerMovie : Plays video with active tracks');
%if isStack == 0 || isStack == 2 %directory with images (allows for saving movie frames)
%    PredTrackerMovie(im_path,vtracks);
%else % .TIFF stack
%    PredTrackerMovie(stack_path,vtracks);
%end