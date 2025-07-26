% PlotTracks.m provides example code to visualize the vtracks
% output of the PredictiveTracker.m function using built-in plotting
% functionalities of the Velocities.m, Vorticities.m, and Divergences.m 
% functions. Additional 
% 
% In order to run, this does require the following:
%   -   Velocities.m : outputs [X, Y, T, U, V, Tr] from vtracks
%                      also plots particle velocity vectors
%   -   Vorticities.m : outputs [Vort, X, Y, T, Tr] from vtracks
%                       also plots particle vorticities
%                       (need to update)
%   -   Divergences.m : outputs [Div, X, Y, T, Tr] from vtracks
%                       also plots particle divergences
%   -   longtracks.m : plots tracks over a given length
%   -   bigread4.m : quickly loads large .TIFF stacks
%   -   TiffTrackViewer.m : viewer for playing images as a movie with tracks overlaid
%   -   PredTrackerMovie.m : plays active tracks over images as a movie --> TO ADD MAYBE

close all
clear all

% Load vtracks structure (assumes you have already run PredictiveTracker.m elsewhere)
%load('vtracks.mat');
load('/Volumes/Sophie/Clam_Flume_Stress_History_Experiments/7_Experiments/Experiment_1/Control/24-09-24_Run1/53_Percent_Pump/output2/8disp_105thresh_vtracks.mat');

% Use built-in plotting functionality of Velocities, Vorticities, and Divergences 
% to visualize particles-----------------------------------------------------------------
framerange=[1 inf];
noisy = 1; %plotting parameter : no plot = 0; plot velocity vectors = 1

disp('Velocities : plots particle velocity vectors');
[u,v,x,y,t,tr]=Velocities(vtracks,framerange,noisy);

disp('Vorticities : plots particle vorticities')
%[vort,~,~,~,~]=Vorticities(vtracks,framerange,noisy) %need to fix old vorticities code

disp('Divergences : plots particle divergences')
%[div,~,~,~,~]=Divergences(vtracks,framerange,noisy) %need to fix old divergences code

% Use longtracks to view tracks over minlength -------------------------------------------
tracks = sortrows([tr x y], 1); %need to sort data trackwise with TR, X, Y
minlength = 20; %minimum track length to plot

fprintf('longtracks : Plots all tracks longer than %d frames.\n', minlength);
[nlong,ntracks,lmean,lrms,lmax]=longtracks(tracks,minlength,noisy);

% Tiff player with functionality for viewing overlaid tracks -----------------------------
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
    %stack_name = 'TrackExampleStack.tiff';
    stack_name = '/Volumes/Sophie/Clam_Flume_Stress_History_Experiments/7_Experiments/Experiment_1/Control/24-09-24_Run1/53_Percent_Pump/Stacked_Bed_Imagery.tiff';
    
    % Parameters for reading .TIFF stacks
    tiff_frame_start = 1; %starting frame of tiff stack
    tiff_num_frames = []; %number of tiff frames to read after starting index
    
    TS = bigread4(stack_name, tiff_frame_start, tiff_num_frames); %load .TIFF stack
end
%%
% Displays viewer with intuitive controls for playing images as a movie and
% overlaying tracks, as well as plotting trajectory metrics
n_ch = 1; %number of channels (1 if greyscale; 3 if color)
minlen = 360;  %minimum track length; if too low, viewer gets very laggy

disp('TiffTrackViewer : TIFF movie viewer and interactive track visualizer');
if isStack == 0 %functionality more limited for directory with images
    tv = TiffTrackViewer('TrackExample', n_ch, vtracks, minlen);
else %otherwise should display .TIFF stack
    tv = TiffTrackViewer(TS,n_ch,vtracks,minlen); 
end

%% Visualization Method 3: Predictive Tracker Movie
% This function plays the same movie as the PredictiveTracker function
%noisy = 1; %noisy = 2 saves indivual frames

%disp('PredTrackerMovie : Plays video with active tracks');
%if isStack == 0 || isStack == 2 %directory with images (allows for saving movie frames)
%    PredTrackerMovie(im_path,vtracks);
%else % .TIFF stack
%    PredTrackerMovie(stack_path,vtracks);
%end