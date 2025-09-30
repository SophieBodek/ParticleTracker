function [atracks,ax,ay,x,y,t,tr]=Accelerations(vtracks,framerange,noisy,order)
% Usage: [atracks,ax,ay,x,y,t,tr]=Accelerations(vtracks,[framerange],[noisy],order)
% Working from the velocity tracks in "vtracks", Accelerations plots the
% acceleration of each particle in frames "framerange". Acceleration values are 
% returned in vectors "ax" and "ay"; corresponding positions, times, and 
% track indices are returned in "x", "y", and "t". If noisy==0, no plot is 
% produced. The input must be a structure of the form produced by 
% PredictiveTracker.m. Specify "framerange" as a two-element vector: 
% [starttime endtime]. Specify method for computing accelerations using
% "order"; order==1 takes first derivative of velocities in vtracks (default),
% while order==2 takes second derivative of positions in vtracks. Function also 
% outputs "Atracks", a structure with the same form as "vtracks" with added
% fields for horizontal and vertical accelerations.

% Based on PredictiveTracker.m and Velocities.m functions
% originally written by Douglas H. Kelley in 2011-12.

framerangedefault = [-inf inf]; % all frames
noisydefault=1;
orderdefault=1;

if nargin<1
    error(['Usage: [atracks,ax,ay,x,y,t,tr] = ' mfilename ...
        '(vtracks,[framerange],[noisy]'])
end
if ~exist('framerange','var') || isempty(framerange)
    framerange=framerangedefault;
elseif numel(framerange)==1
    framerange=framerange*[1 1];
end
if ~exist('noisy','var') || isempty(noisy)
    noisy=noisydefault;
end
if ~exist('order','var') || isempty(order)
    order=orderdefault;
end

fn=fieldnames(vtracks);
if ~any(strcmp(fn,'X')) || ~any(strcmp(fn,'Y')) || ~any(strcmp(fn,'T')) ...
        || ~any(strcmp(fn,'U')) || ~any(strcmp(fn,'V')) || isempty(vtracks)
    error('Sorry, the input does not appear to contain tracks.')
end
if any(strcmp(fn,'Theta')) % if there is a field named Theta
    minarea = 2; %particles must have had a minarea =/ 1
else 
    minarea = 1;
end

% Re-assemble vtracks structure, except with accelerations calculated ------------------
filterwidth = 1; fitwidth = 3; % parameters for the differentiation kernel
% Prune lengths --> recall derivatives will shorten tracks
vtracks = vtracks([vtracks.len] >= (2*fitwidth+1));
ntracks = numel(vtracks);

% convolution kernel - first derivative
Av = 1.0/(0.5*filterwidth^2 * ...
    (sqrt(pi)*filterwidth*erf(fitwidth/filterwidth) - ...
    2*fitwidth*exp(-fitwidth^2/filterwidth^2)));
vkernel = -fitwidth:fitwidth;
vkernel = Av.*vkernel.*exp(-vkernel.^2./filterwidth^2); 

% Normalization factor for second derivative
Aa = 2.0/(filterwidth^4 * ...
    (sqrt(pi)*filterwidth*erf(fitwidth/fitwidth) - ...
    2*fitwidth*exp(-fitwidth^2/filterwidth^2) + ...
    sqrt(pi)*filterwidth^3*erf(fitwidth/filterwidth)/2));

akernel = -fitwidth:fitwidth;
akernel = Aa .* (2*akernel.^2/filterwidth^2 - 1) .* exp(-akernel.^2./filterwidth^2);

 % loop over tracks
if minarea==1
    atracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[],'Ax',[],'Ay',[]),ntracks,1);
else
    atracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[], ...
        'V',[],'Ax',[],'Ay',[],'Theta',[]),ntracks,1);
end

for ii = 1:ntracks
    if order == 1
        % Compute accelerations (first derivative of velocity)
        ax = -conv(vtracks(ii).U, vkernel, 'valid'); 
        ay = -conv(vtracks(ii).V, vkernel, 'valid');  
    elseif order == 2
        % Compute accelerations (second derivative of position)
        ax = conv(vtracks(ii).X,akernel,'valid');  
        ay = conv(vtracks(ii).Y,akernel,'valid');
    end
    if minarea==1
        atracks(ii) = struct('len',vtracks(ii).len - 2*fitwidth, ...
            'X',vtracks(ii).X(fitwidth+1:end-fitwidth), ...
            'Y',vtracks(ii).Y(fitwidth+1:end-fitwidth), ...
            'T',vtracks(ii).T(fitwidth+1:end-fitwidth), ...
            'U',vtracks(ii).U(fitwidth+1:end-fitwidth), ...
            'V',vtracks(ii).V(fitwidth+1:end-fitwidth), ...
            'Ax',ax, ...
            'Ay',ay);
    else
        atracks(ii) = struct('len',vtracks(ii).len - 2*fitwidth, ...
            'X',vtracks(ii).X(fitwidth+1:end-fitwidth), ...
            'Y',vtracks(ii).Y(fitwidth+1:end-fitwidth), ...
            'T',vtracks(ii).T(fitwidth+1:end-fitwidth), ...
            'U',vtracks(ii).U(fitwidth+1:end-fitwidth), ...
            'V',vtracks(ii).V(fitwidth+1:end-fitwidth), ...
            'Ax',ax, ...
            'Ay',ay, ...
            'Theta',vtracks(ii).Theta(fitwidth+1:end-fitwidth));
    end
end


 % -----------------------------------------------------
Nall=sum([atracks.len]);
ax=NaN(Nall,1); % initialize arrays
ay=NaN(Nall,1);
x=NaN(Nall,1);
y=NaN(Nall,1);
t=NaN(Nall,1);
tr=NaN(Nall,1);
pos=1; % marker to current position in arrays
for ii=1:numel(atracks) % assemble arrays from atracks structure
    ind = (atracks(ii).T>=min(framerange)) & (atracks(ii).T<=max(framerange));
    if any(ind)
        Np=sum(ind);
        ax(pos:pos+Np-1)=atracks(ii).Ax(ind(:));
        ay(pos:pos+Np-1)=atracks(ii).Ay(ind(:));
        x(pos:pos+Np-1)=atracks(ii).X(ind(:));
        y(pos:pos+Np-1)=atracks(ii).Y(ind(:));
        t(pos:pos+Np-1)=atracks(ii).T(ind(:));
        tr(pos:pos+Np-1)=ii;
        pos=pos+Np;
    end
end
ind=~isnan(ax); % remove unwanted frames
ax=ax(ind);
ay=ay(ind);
x=x(ind);
y=y(ind);
t=t(ind);
tr=tr(ind);
[t,ind]=sort(t(:)); % sort by time
ax=ax(ind);
ay=ay(ind);
x=x(ind);
y=y(ind);
tr=tr(ind);

% If requested, plot the particles
if noisy
    figure;
    
    % Calculate both magnitude and signed components
    accel_magnitude = sqrt(ax.^2 + ay.^2);
    accel_signed = sign(ax + ay) .* accel_magnitude; %sign based on dominant direction
    
    % Plot
    scatter(x, y, 60, accel_signed, 'filled');
    hold on;
    quiver(x, y, ax, ay, 0.5, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.3);
    colormap(flipud(cool)); % blue-cyan-magenta-red
    
    % Center colorbar at zero
    max_val = max(abs(accel_signed));
    if max_val > 0
        clim([-max_val, max_val]);
    end

    cb = colorbar;
    ylabel(cb, 'Sign of corresponds to dominant direction');
    set(gca,'dataaspectratio',[1 1 1]);
    axis tight
    if framerange(1)==framerange(2)
        title(['vorticity at ' num2str(numel(x)) ' particles in frame ' ...
            num2str(framerange(1))]);
    else
        if isinf(framerange(1))
            framerange(1)=t(1);
        end
        if isinf(framerange(2))
            framerange(2)=t(end);
        end
        title(['accelerations at ' num2str(numel(x)) ' particles in frames ' ...
            num2str(min(framerange)) ' to ' num2str(max(framerange)) ]);
    end
end

