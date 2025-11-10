function [ktracks,k,x,y,t,tr]=Curvature(vtracks,framerange,noisy)
% Usage: [ktracks,k,x,y,t,tr]=Curvature(vtracks,[framerange],[noisy])
% Working from the velocity and acceleration tracks in "vtracks", Curvature calculates
% curvatures of the trajectories within the frames in "framerange". Curvature is 
% returned in an added field "K" in the ktracks struct. Curvature values are 
% also returned in the vector "k"; corresponding positions, times, and track 
% indices are returned in "x", "y", "t", and "tr". Specify "framerange" as a 
% two-element vector: [starttime endtime]. If noisy==0, no plot is produced. 
%
% The input must be a structure of the form produced by PredictiveTracker.m, and 
% requires fields for U, V, Ax, and Ay.

% Based on PredictiveTracker.m and Velocities.m functions
% originally written by Douglas H. Kelley in 2011-12.

framerangedefault = [-inf inf]; % all frames
noisydefault=1; 
thresh_speed = 0.01; % filter out low speeds

if nargin<1
    error(['Usage: [ktracks,k,x,y,t,tr] = ' mfilename ...
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

fn=fieldnames(vtracks); % needs U, V, Ax, and Ay
if ~any(strcmp(fn,'X')) || ~any(strcmp(fn,'Y')) || ~any(strcmp(fn,'T')) ...
        || ~any(strcmp(fn,'U')) || ~any(strcmp(fn,'V')) ...
        || ~any(strcmp(fn,'Ax')) || ~any(strcmp(fn,'Ay')) || isempty(vtracks)
    error('Sorry, the input does not appear to contain tracks.')
end
if any(strcmp(fn,'Theta')) % if there is a field named Theta
    minarea = 2; %particles must have had a minarea =/ 1
else 
    minarea = 1;
end

% initialize ktracks with additional curvature field
ntracks = length(vtracks);
if minarea==1
katracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[],'V',[],'Ax',[],'Ay',[],'K',[]),ntracks,1);
else
    ktracks = repmat(struct('len',[],'X',[],'Y',[],'T',[],'U',[], ...
            'V',[],'Ax',[],'Ay',[],'K',[],'Theta',[]),ntracks,1);
end

% loop over tracks
for ii = 1:ntracks
    % Calculate curvature
    U = vtracks(ii).U; 
    V = vtracks(ii).V;
    Ax = vtracks(ii).Ax; 
    Ay = vtracks(ii).Ay;

    speed = sqrt(U.^2 + V.^2);
    kappa = abs(U.*Ay - V.*Ax) ./ (speed.^3 + eps);
    kappa(speed < thresh_speed) = 0;  % Mask low speeds -- what's a good value for this!?!?

    % Build and add curvature to atracks struct
    if minarea==1
        ktracks(ii) = struct('len',length(kappa), ... 
            'X',vtracks(ii).X, ...
            'Y',vtracks(ii).Y, ...
            'T',vtracks(ii).T, ...
            'U',vtracks(ii).U, ...
            'V',vtracks(ii).V, ...
            'Ax',vtracks(ii).Ax, ...
            'Ay',vtracks(ii).Ay,...
            'K', kappa);
    else
        ktracks(ii) = struct('len',length(kappa), ... 
            'X',vtracks(ii).X, ...
            'Y',vtracks(ii).Y, ...
            'T',vtracks(ii).T, ...
            'U',vtracks(ii).U, ...
            'V',vtracks(ii).V, ...
            'Ax',vtracks(ii).Ax, ...
            'Ay',vtracks(ii).Ay,...
            'Theta',vtracks(ii).Theta,...
            'K', kappa);
    end
end


% Assemble output vectors: ax, ay, x, y, t, tr, [curv] ----------------------
Nall=sum([ktracks.len]);
ax=NaN(Nall,1); % initialize arrays
ay=NaN(Nall,1);
x=NaN(Nall,1);
y=NaN(Nall,1);
t=NaN(Nall,1);
tr=NaN(Nall,1);
k=NaN(Nall,1);

pos=1; % marker to current position in arrays
for ii=1:numel(ktracks) % assemble arrays from atracks structure
    ind = (ktracks(ii).T>=min(framerange)) & (ktracks(ii).T<=max(framerange));
    if any(ind)
        Np=sum(ind);
        %ax(pos:pos+Np-1)=ktracks(ii).Ax(ind(:));
        %ay(pos:pos+Np-1)=ktracks(ii).Ay(ind(:));
        %u(pos:pos+Np-1)=ktracks(ii).U(ind(:));
        %v(pos:pos+Np-1)=ktracks(ii).V(ind(:));
        x(pos:pos+Np-1)=ktracks(ii).X(ind(:));
        y(pos:pos+Np-1)=ktracks(ii).Y(ind(:));
        t(pos:pos+Np-1)=ktracks(ii).T(ind(:));
        tr(pos:pos+Np-1)=ii;
        k(pos:pos+Np-1)=ktracks(ii).K(ind(:));
        pos=pos+Np;
    end
end

ind=~isnan(k); % remove unwanted frames
ax=ax(ind);
ay=ay(ind);
%u=u(ind);
%v=v(ind);
x=x(ind);
y=y(ind);
t=t(ind);
tr=tr(ind);
k=k(ind);

[t,ind]=sort(t(:)); % sort by time
ax=ax(ind);
ay=ay(ind);
%u=u(ind);
%v=v(ind);
x=x(ind);
y=y(ind);
tr=tr(ind);
k=k(ind);

% If requested, plot the particles
if noisy
    figure;
    
    % Plot curvature with colormap
    scatter(x, y, 30, k, 'filled');
     hold on;
    quiver(x, y, ax, ay, 0.5, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.3);
        
    colormap(cool);
    cb = colorbar;
    ylabel(cb, 'Curvature (1/[Length Unit])', 'FontSize', 12);
        
    set(gca, 'DataAspectRatio', [1 1 1]);
    axis tight;
    xlabel('X [units]', 'FontSize', 12);
    ylabel('Y [units]', 'FontSize', 12);
        
    % Title based on frame range
    if framerange(1)==framerange(2)
        title(['curvature of ' num2str(numel(x)) ' particles in frame ' ...
               num2str(framerange(1))], 'FontSize', 14);
    else
        if isinf(framerange(1))
            framerange(1)=t(1);
        end
        if isinf(framerange(2))
            framerange(2)=t(end);
        end
        title(['curvature of ' num2str(numel(x)) ' particles in frames ' ...
               num2str(min(framerange)) ' to ' num2str(max(framerange))]);
        end
end

end
