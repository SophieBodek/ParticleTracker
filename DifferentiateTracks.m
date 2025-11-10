function vtracks = DifferentiateTracks(tracks,filterwidth,fitwidth,smoothed)
% Usage : [vtracks] = DifferentiateTracks(tracks,filterwidth,fitwidth,[smoothed])
%
% Calculates particle velocities and accelerations from trajectories using
% convolution with a smoothing and differentiating kernel. Method is 
% described in Mordant et al. (2004) and the calculation is fully outlined
% in Ouelette (2006; PhD thesis). The DifferentiateTracks function was 
% originally written in C++ by N.T. Ouelette in 8/1/2008. 
% 
% Inputs:
%   tracks      - struct array output from PredictiveTracker with 
%                 fields: len, X, Y, T, Theta (optional)
%   filterwidth - kernel width (recommend filterwidth = fitwidth/1.5)
%   fitwidth    - width of sliding window (full size 2*fitwidth + 1) 
%   smoothed    - option to calculate smoothed trajectories (smoothed==1);
%                 otherwise only outputs velocities and accelerations
%                 (smoothed==0; default)
%
% Outputs:
%   vtracks     - struct array with fields: len, X, Y, T, U, V, Ax, Ay, 
%                 Xfilt, Yfilt (if smoothed==1), Theta (if originally included)

% -=- SETUP -=------------------------------------------------------------------------
if nargin<1
    error(['Usage: vtracks = ' mfilename '(tracks,filterwidth,fitwidthl,[smoothed]'])
end
if ~exist('smoothed','var') || isempty(smoothed)
    smoothed=0; % default is not include smoothed trajectories
end

fn=fieldnames(tracks);
if ~any(strcmp(fn,'X')) || ~any(strcmp(fn,'Y')) || ~any(strcmp(fn,'T')) ...
        || ~any(strcmp(fn,'len')) || isempty(tracks)
    error('Sorry, the input does not appear to contain tracks.')
end

% Check whether tracks struct has Theta field
if isfield(tracks, 'Theta')
    minarea = 2;
else
    minarea = 1; 
end

% -=- Prune tracks that are too short -=------------------------------------------
disp('Pruning...');
tracks = tracks([tracks.len] > (2*fitwidth+1));
ntracks = numel(tracks);
%meanlength = mean([tracks.len]);
%rmslength = sqrt(mean([tracks.len].^2));

% -=- Kernels -=------------------------------------------------------------------------
t_vals = -fitwidth:fitwidth;

% define the convolution kernel for smoothed positions
Bp = exp(-fitwidth^2/filterwidth^2) / (2*fitwidth*exp(-fitwidth^2/filterwidth^2) - ...
     filterwidth*sqrt(pi)*erf(fitwidth/filterwidth));
Ap = (1 - 2*fitwidth*Bp) / (filterwidth*sqrt(pi)*erf(fitwidth/filterwidth));
pkernel = Ap .* exp(-t_vals.^2 / filterwidth^2) + Bp;

% define the convolution kernel for velocities
Av = 1.0/(0.5*filterwidth^2 * ...
    (sqrt(pi)*filterwidth*erf(fitwidth/filterwidth) - ...
    2*fitwidth*exp(-fitwidth^2/filterwidth^2)));
vkernel = Av .* t_vals .* exp(-t_vals.^2 / filterwidth^2);

% define the convolution kernel for accelerations
[Aa, Ba] = CalcConst(fitwidth, filterwidth);
akernel = Aa .* (2 .* t_vals.^2 / filterwidth^2 - 1) .* exp(-t_vals.^2 / filterwidth^2) + Ba;

% -=- Differetiate -=-------------------------------------------------------
disp('Smoothing & differentiating...');

% Set up data struct
if minarea==1 && smoothed==1
    vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[],'Ax',[],'Ay',[],...
                     'Xfilt',[],'Yfilt',[]),ntracks,1);
elseif minarea==2 && smoothed==1
    vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[], ...
                     'V',[],'Ax',[],'Ay',[],'Xfilt',[],'Yfilt',[],'Theta',[]),ntracks,1);
elseif minarea==1 && smoothed==0
    vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[],'Ax',[],'Ay',[]),ntracks,1);
else %minarea==2 && smoothed==0
    vtracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[],'Ax',[],'Ay',[],'Theta',[]),ntracks,1);
end

% Loop over tracks
for ii = 1:ntracks
    tracklength = tracks(ii).len;
        
    % Check if track is too short (pruning should eliminate)
    if tracklength < (2*fitwidth + 1)
        continue;
    end

    % Get time step
    dT = tracks(ii).T(2) - tracks(ii).T(1);
        
    % Number of output points
    npts = tracklength - 2*fitwidth - 1;
        
    % Preallocate output arrays
    X_out = zeros(npts, 1)';
    Y_out = zeros(npts, 1)';
    T_out = zeros(npts, 1)';
    U_out = zeros(npts, 1)';
    V_out = zeros(npts, 1)';
    Ax_out = zeros(npts, 1)';
    Ay_out = zeros(npts, 1)';
    if smoothed ==1 %include smoothed tracks
        Xfilt_out = zeros(npts,1)';
        Yfilt_out = zeros(npts,1)';
    end
    if minarea == 2 %has angle/Theta
        Theta_out = zeros(npts, 1)';
    end
        
    % Loop over the track with sliding window of 2*fitwidth+1
    for jj = (fitwidth+1):(tracklength - fitwidth-1)
        idx = jj - fitwidth;
            
        % Extract window of positions
        window_idx = (jj-fitwidth):(jj+fitwidth);
        valsx = tracks(ii).X(window_idx);
        valsy = tracks(ii).Y(window_idx);
            
        % Store central position and time
        X_out(idx) = tracks(ii).X(jj);
        Y_out(idx) = tracks(ii).Y(jj);
        T_out(idx) = tracks(ii).T(jj);
        if minarea==2
            Theta_out(idx) = tracks(ii).Theta(jj);
        end

        % Compute smoothed position (convolution with position kernel)
        if smoothed ==1
            Xfilt_out(idx) = Convolve(valsx, pkernel, fitwidth);
            Yfilt_out(idx) = Convolve(valsy, pkernel, fitwidth);
        end
            
        % Compute velocity (convolution with velocity kernel)
        U_out(idx) = Convolve(valsx, vkernel, fitwidth) / dT;%-conv(valsx,vkernel,'valid');
        V_out(idx) = Convolve(valsy, vkernel, fitwidth) / dT;%-conv(valsy,vkernel,'valid');
            
        % Compute acceleration (convolution with acceleration kernel)
        Ax_out(idx) = Convolve(valsx, akernel, fitwidth) / (dT * dT);
        Ay_out(idx) = Convolve(valsy, akernel, fitwidth) / (dT * dT);

    end
        
    % Store results
    if minarea ==1 && smoothed==1
        vtracks(ii) = struct('len', npts, 'X', X_out, 'Y', Y_out, ...
                            'T', T_out, 'U', U_out, 'V', V_out, ...
                            'Ax', Ax_out, 'Ay', Ay_out, ...
                            'Xfilt', Xfilt_out, 'Yfilt', Yfilt_out);
    elseif minarea==2 && smoothed==1
        vtracks(ii) = struct('len', npts, 'X', X_out, 'Y', Y_out, ...
                            'T', T_out, 'U', U_out, 'V', V_out, ...
                            'Ax', Ax_out, 'Ay', Ay_out, ...
                            'Xfilt', Xfilt_out, 'Yfilt', Yfilt_out, 'Theta', Theta_out);
    elseif minarea==1 && smoothed==0
        vtracks(ii) = struct('len', npts, 'X', X_out, 'Y', Y_out, ...
                             'T', T_out, 'U', U_out, 'V', V_out, ...
                            'Ax', Ax_out, 'Ay', Ay_out);
    else %minarea==2 && smoothed==0
        vtracks(ii) = struct('len', npts, 'X', X_out, 'Y', Y_out, ...
                            'T', T_out, 'U', U_out, 'V', V_out, ...
                            'Ax', Ax_out, 'Ay', Ay_out, 'Theta', Theta_out);
    end
        
    % Progress update
    if mod(ii, 100) == 0
        fprintf('Processed %d of %d tracks.\n', ii, ntracks);
    end
end
    
disp('Done.');
end

% -=- HELPER FUNCTIONS -=-----------------------------------------------
function [Aa, Ba] = CalcConst(L, w)
    % Calculate constants for acceleration kernel
    sumA = 0;
    sumB = 0;
    
    for i = 1:(2*L - 1)
        t = -L + i;
        
        tmpA = (0.5 * t^2 - L^2/6.0) * (2.0 * t^2 / w^2 - 1.0) * exp(-t^2 / w^2);
        tmpB = (2.0 * t^2 / w^2 - 1.0) * exp(-t^2 / w^2);
        
        if mod(i, 2) == 1  % odd
            sumA = sumA + tmpA * 4.0;
            sumB = sumB + tmpB * 4.0;
        else  % even
            sumA = sumA + tmpA * 2.0;
            sumB = sumB + tmpB * 2.0;
        end
    end
    
    % Add endpoints
    t = -L;
    sumA = sumA + (0.5 * t^2 - L^2/6.0) * (2.0 * t^2 / w^2 - 1.0) * exp(-t^2 / w^2);
    sumB = sumB + (2.0 * t^2 / w^2 - 1.0) * exp(-t^2 / w^2);
    
    t = L;
    sumA = sumA + (0.5 * t^2 - L^2/6.0) * (2.0 * t^2 / w^2 - 1.0) * exp(-t^2 / w^2);
    sumB = sumB + (2.0 * t^2 / w^2 - 1.0) * exp(-t^2 / w^2);
    
    sumA = sumA / 3.0;
    sumB = sumB / 3.0;
    
    Aa = 1.0 / sumA;
    Ba = (-1.0 / (2.0 * L)) * Aa * sumB;
end

function sum_val = Convolve(vals, kvals, L)
    sum_val = 0;
    
    for i = 2:(2*L)  
        if mod(i, 2) == 0 %odd
            sum_val = sum_val + vals(i) * kvals(i) * 4.0;
        else 
            sum_val = sum_val + vals(i) * kvals(i) * 2.0;
        end
    end
    
    % Add the first and last elements
    sum_val = sum_val + vals(1) * kvals(1);
    sum_val = sum_val + vals(2*L + 1) * kvals(2*L + 1);
    
    sum_val = sum_val / 3.0;
end