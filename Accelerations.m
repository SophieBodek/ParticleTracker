function [ax,ay,x,y,t,tr]=Accelerations(vtracks,framerange,noisy)
% Usage: [ax,ay,x,y,t,tr]=Accelerations(vtracks,[framerange],[noisy])
% Working from the acceleration tracks in "vtracks", Accelerations plots
% the acceleration of each particle in frames "framerange". Acceleration values 
% are returned in vectors "ax" and "ay"; corresponding positions, times, and 
% track indices are returned in "x", "y", and "t". Specify "framerange" as a 
% two-element vector: [starttime endtime]. If noisy==0, no plot is produced. 
%
% The input must be a structure of the form produced by DifferentiateTracks.m
% (or PredictiveTracker.m) with at least fields "X", "Y", "T", "Ax", and "Ay".
%
% Based on PredictiveTracker.m and Velocities.m functions
% originally written by Douglas H. Kelley in 2011-12.

framerangedefault = [-inf inf]; % all frames
noisydefault=1; 

if nargin<1
    error(['Usage: [ax,ay,x,y,t,tr] = ' mfilename ...
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

fn=fieldnames(vtracks); % needs position, time, and accelerations
if ~any(strcmp(fn,'X')) || ~any(strcmp(fn,'Y')) || ~any(strcmp(fn,'T')) ...
        || ~any(strcmp(fn,'Ax')) || ~any(strcmp(fn,'Ay')) || isempty(vtracks)
    error('Sorry, the input does not appear to contain tracks.')
end

% Assemble output vectors: ax, ay, x, y, t, tr -------------------------
Nall=sum([vtracks.len]);
ax=NaN(Nall,1); % initialize arrays
ay=NaN(Nall,1);
x=NaN(Nall,1);
y=NaN(Nall,1);
t=NaN(Nall,1);
tr=NaN(Nall,1);

pos=1; % marker to current position in arrays
for ii=1:numel(vtracks) % assemble arrays from atracks structure
    ind = (vtracks(ii).T>=min(framerange)) & (vtracks(ii).T<=max(framerange));
    if any(ind)
        Np=sum(ind);
        ax(pos:pos+Np-1)=vtracks(ii).Ax(ind(:));
        ay(pos:pos+Np-1)=vtracks(ii).Ay(ind(:));
        %u(pos:pos+Np-1)=vtracks(ii).U(ind(:));
        %v(pos:pos+Np-1)=vtracks(ii).V(ind(:));
        x(pos:pos+Np-1)=vtracks(ii).X(ind(:));
        y(pos:pos+Np-1)=vtracks(ii).Y(ind(:));
        t(pos:pos+Np-1)=vtracks(ii).T(ind(:));
        tr(pos:pos+Np-1)=ii;
        pos=pos+Np;
    end
end

ind=~isnan(ax); % remove unwanted frames
ax=ax(ind);
ay=ay(ind);
%u=u(ind);
%v=v(ind);
x=x(ind);
y=y(ind);
t=t(ind);
tr=tr(ind);
[t,ind]=sort(t(:)); % sort by time
ax=ax(ind);
ay=ay(ind);
%u=u(ind);
%v=v(ind);
x=x(ind);
y=y(ind);
tr=tr(ind);

% If requested, plot the particles
if noisy
    figure; %plot accelerations with colormap
    accel_magnitude = sqrt(ax.^2 + ay.^2);
    accel_signed = sign(ax + ay) .* accel_magnitude; %sign based on dominant direction
        
    scatter(x, y, 40, accel_signed, 'filled');
    hold on;
    quiver(x, y, ax, ay, 0.5, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.3);
    colormap(flipud(turbo)); 
        
    % Center colorbar at zero
    max_val = max(abs(accel_signed));
    if max_val > 0
        clim([-max_val, max_val]);
    end

    cb = colorbar;
    ylabel(cb, 'Acceleration magnitude (signed) [L/T^2]', 'FontSize', 12);
    set(gca, 'DataAspectRatio', [1 1 1]);
    axis tight;
    xlabel('X [units]', 'FontSize', 12);
    ylabel('Y [units]', 'FontSize', 12);
        
     if framerange(1)==framerange(2)
        title(['accelerations at ' num2str(numel(x)) ' particles in frame ' ...
            num2str(framerange(1))], 'FontSize', 14);
     else
        if isinf(framerange(1))
            framerange(1)=t(1);
        end
        if isinf(framerange(2))
            framerange(2)=t(end);
        end
        title(['accelerations at ' num2str(numel(x)) ' particles in frames ' ...
            num2str(min(framerange)) ' to ' num2str(max(framerange))]);
    end
end

end

