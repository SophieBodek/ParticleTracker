% TrackExample.m, written 13 February 2012 by Doug Kelley, walks the user
% through the use of the particle tracking function PredictiveTracker.m 
% and the associated functions BackgroundImage.m and Velocities.m. This file 
% can be downloaded from https://web.stanford.edu/~nto/software_tracking.shtml.

clear all

fprintf('\n');
disp('Welcome to particle tracking! This script, "TrackExample.m" will walk ')
disp('you through a simple example of tracking particles using BackgroundImage.m, ')
disp('PredictiveTracker.m, and Velocities.m. All have been developed by ')
disp('Nicholas T. Ouellette (Stanford University) and Douglas H. Kelley (University of ')
disp('Rochester) and are available from https://web.stanford.edu/~nto/software_tracking.shtml. ')
input('Press enter to continue. ')

fprintf('\n');
disp('PredictiveTracker is written to track particles in movies saved as image ')
disp('sequences. That is, each frame is saved as a separate image file, and their ')
disp('order is given by their filenames. In the same directory as this script you ')
disp('should find 100 tiff images saved in the TrackExample folder, with names ranging ')
disp('from "0004.tif" to "0103.tif". Each image shows fluorescent tracer particles in ')
disp('a fluid experiment conducted at Yale University in the Ouellette lab. You can view ')
disp('the frames with any image software you like. Press enter to view the first frame ') 
disp('by running these Matlab commands: ')
disp('   figure; fr1=imread(''TrackExample/0004.tif''); imagesc(fr1); colormap gray; ')
input('   title(''0004.tif''); ')
figure; fr1=imread('TrackExample/0004.tif'); imagesc(fr1); colormap gray; title('0004.tif');

fprintf('\n');
disp('Now, let''s make a background image. We don''t want to track anything that ')
disp('doesn''t move, so we''ll calculate the ensemble mean image of all the ')
disp('frames of our movie (0*.tif). Moving objects will tend to get smeared ')
disp('out and disappear, while stationary objects will remain distinct. (Later ')
disp('PredictiveTracker will subtract the background image from each frame before tracking.) ')
disp('We''ll save the background images as "background.tif" in the TrackExample folder. ')
disp('Press enter to run this command: ')
input('   bg = BackgroundImage(''TrackExample/0*.tif'',''TrackExample/background.tif''); ')
bg = BackgroundImage('TrackExample/0*.tif','TrackExample/background.tif');

fprintf('\n');
disp('Let''s see what the background image looks like. ')
disp('Press enter to run these commands: ')
input('   figure; imagesc(bg); colormap gray; title(''background image''); ');
figure; imagesc(bg); colormap gray; title('background image');

fprintf('\n');
disp('For this example we''re using only 100 frames, so the background image is')
disp('a bit streaky. With more frames, moving objects will disappear much more ')
disp('fully. Now we''ll use PredictiveTracker to do the tracking. Entering the ')
disp('name of the function without any inputs results in a syntax prompt. Press ')
disp('enter to run this command:')
input('   PredictiveTracker ')
try
    PredictiveTracker
catch err
    fprintf('\n');
    disp(['??? Error using ==> PredictiveTracker at ' num2str(err.stack.line)]) % mimic usual output
    beep
    disp(err.message)
end

fprintf('\n');
disp('For this example, we''ll choose a threshold of 10 -- that is, to count as a ')
disp('potential particle, a pixel must be brighter than its neighbors and at least ')
disp('10 bits brighter than the background. And we''ll choose max_disp = 8 -- that ')
disp('is, to continue a track, a potential particle must lie within 8 pixels of the ')
disp('predicted particle location. We''ll use the background image created above ')
disp('(''background.tif''). We''ll choose minarea = 1 -- that is, we expect ')
disp('particles to be single pixels. In this movie particles are bright, so we ')
disp('don''t want to invert (invert = 0), and we''ll ask to see a movie showing ')
disp('the results (noisy = 1). The outputs are "vtracks" (a struct array containing ')
disp('the tracks), "ntr" (the number of tracks), "lm" (the mean track length), and ')
disp('"lrms" (the root-mean-square track length). PredictiveTracker will report its ')
disp('progress as it goes. Press enter to run this command:') 
disp('   [vtracks,ntr,lm,lrms] = ...')
input('       PredictiveTracker(''TrackExample/0*.tif'',10,8,''TrackExample/background.tif'',1,0,1); ');
[vtracks,ntr,lm,lrms] = ...
    PredictiveTracker('TrackExample/0*.tif',10,8,'TrackExample/background.tif',1,0,1);
fprintf('\n');

disp('A figure appears as requested, and the movie plays back with particle tracks ')
disp('overlaid in color. It''s a useful way to check tracking parameters. If too ')
disp('few particles are identified, consider lowering the threshold. If too many ')
disp('particles are identified and the code runs very slowly, consider raising the ')
disp('threshold. If tracks are too short, consider raising max_disp. If tracks are ')
disp('noisy or make erroneous connections, consider lowering max_disp. To track ')
disp('without seeing the movie, use noisy = 0. ')

fprintf('\n');
disp('PredictiveTracker outputs tracks as a struct array with one element per ')
input('track. Press enter to see the form of the "vtracks" struct array. ')
vtracks

fprintf('\n');
disp('Each element contains a scalar giving the track length, as well as vectors ')
disp('giving the positions, velocities, and times of the track. All units are in ')
input('terms of pixels and frames. Press enter to see the form of the first track. ');
vtracks(1)

fprintf('\n');
disp('Now let''s plot the results another way, using the Velocities function. Press ')
disp('enter to see a syntax prompt by running Velocities without any inputs: ')
input('   Velocities ')
try
    Velocities
catch err
    fprintf('\n');
    disp(['??? Error using ==> Velocities at ' num2str(err.stack.line)]) % mimic usual output
    beep
    disp(err.message)
end

fprintf('\n');
disp('We''ll ask for data from all frames ( framerange = [0 inf] ), and we do want a ')
disp('plot (noisy = 1). Press enter to run this command: ')
input('   [u,v,x,y,t]=Velocities(vtracks,[0 inf],1); ')
[u,v,x,y,t]=Velocities(vtracks,[0 inf],1);

fprintf('\n');
disp('If you zoom in, you''ll see that each particle is represented by a ')
disp('velocity arrow. Besides producing a plot, Velocities returns ')
disp('positions, velocities, and times as vectors sorted time-wise (not track- ')
disp('wise) which is often convenient. Press enter to calculate the root- ')
disp('mean-square velocity with this command: ')
input('   Urms = mean(sqrt( u.^2+v.^2 )) ');
Urms = mean(sqrt( u.^2+v.^2 ))

fprintf('\n');
disp('That''s the end of our example. The tracks you created are still in memory, ')
disp('so feel free to experiment and analyze the data more. The current version of this ')
disp('software can be downloaded from https://web.stanford.edu/~nto/software_tracking.shtml. ')
disp('This software is freely available; we ask that it be credited when modified ')
disp('or passed along, and that publications resulting from the use of this software ')
disp('cite D. H. Kelley and N. T. Ouellette, "Using particle tracking to measure ')
disp('flow instabilities in an undergraduate laboratory experiment," Am. J. Phys. 79: ')
disp('267-273 (2011). Enjoy! ')
fprintf('\n');

