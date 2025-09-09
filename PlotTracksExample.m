% PlotTracksExample.m provides example code to visualize the vtracks data
% structure of the PredictiveTracker.m function using built-in plotting
% functionalities of the Velocities.m, Vorticities.m, and Divergences.m 
% functions. 
% 
% In order to run, this does require the following:
%   -   Velocities.m : outputs [X, Y, T, U, V, Tr] from vtracks
%                      also plots particle velocity vectors
%   -   Accelerations.m : outputs [Atracks, Ax, Ay, X, Y, T, Tr] from vtracks
%                         also plots particle accelerations
%   -   Vorticities.m : outputs [Vort, X, Y, T, Tr] from vtracks
%                       also plots particle vorticities
%                       (requires Partial Differential Equation toolbox)
%   -   Divergences.m : outputs [Div, X, Y, T, Tr] from vtracks
%                       also plots particle divergences
%                        (requires Partial Differential Equation toolbox)
%   -   longtracks.m : plots tracks over a given length
%   -   bigread4.m : quickly loads large .TIFF stacks
%   -   TiffTrackViewer.m : viewer for playing images as a movie with tracks overlaid
%   -   PredTrackerMovie.m : plays active tracks over images as a movie --> TO ADD MAYBE

%close all
%clear all

% Run background and predictivetracker
%bg = BackgroundImage('TrackExample/0*.tif','TrackExample/background.tif');
%[vtracks_acc,ntracks,lmean,lrms] = ...
%    PredictiveTracker('TrackExample/0*.tif',10,8,'TrackExample/background.tif',1,0,0); 

% OR load vtracks structure (assumes you have already run PredictiveTracker.m elsewhere)
load('vtracks.mat');

% Use built-in plotting functionality of Velocities, Accelerations, Vorticities, and Divergences 
% to visualize particles-----------------------------------------------------------------
framerange=[1 inf];
noisy = 1; %plotting parameter : no plot = 0; plot velocity vectors = 1

disp('Velocities : plots particle velocity vectors');
[u,v,x,y,t,tr]=Velocities(vtracks,framerange,noisy);

disp('Accelerations : plots particle accelerations');
[atracks,ax,ay,~,~,~,~]=Accelerations(vtracks,framerange,noisy,1);
% Note: x,y,t,tr vectors outputted by Accelerations will be shorter than
% those outputted by Velocities due to convolution and smoothing calculation

disp('Vorticities : plots particle vorticities')
[vort,~,~,~,~]=Vorticities(vtracks,framerange,noisy); 

disp('Divergences : plots particle divergences')
[div,~,~,~,~]=Divergences(vtracks,framerange,noisy); 

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
    stack_name = 'TrackExampleStack.tiff';

    % Parameters for reading .TIFF stacks
    tiff_frame_start = 1; %starting frame of tiff stack
    tiff_num_frames = []; %number of tiff frames to read after starting index

    TS = bigread4(stack_name, tiff_frame_start, tiff_num_frames); %load .TIFF stack
end

% Displays viewer with intuitive controls for playing images as a movie and
% overlaying tracks, as well as plotting trajectory metrics
n_ch = 1; %number of channels (1 if greyscale; 3 if color)
minlen = 20;  %minimum track length; if too low, viewer gets very laggy

disp('TiffTrackViewer : TIFF movie viewer and interactive track visualizer');
if isStack == 0 %functionality more limited for directory with images
    tv = TiffTrackViewer('TrackExample', n_ch, atracks, minlen);
else %otherwise should display .TIFF stack
    tv = TiffTrackViewer(TS,n_ch,atracks,minlen); 
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