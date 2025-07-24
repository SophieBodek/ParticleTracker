function [nlong,ntracks,lmean,lrms,lmax]=longtracks(filename,minlength,noisy,readwhole)
% Usage: [nlong,ntracks,lmean,lrms,lmax]=longtracks(filename,[minlength],[noisy],[readwhole])
% Working from the tracks in filename, longtracks plots all tracks
% longer than minlength. It also returns simple statistics about track
% lengths. If noisy=0, no plot is produced. If readwhole==0, a slower, 
% less memory-intensive plotting algorithm is used (will eliminate some
% out-of-memory errors). The input file must be sorted trackwise and have 
% trackindex, x, and y as its first three columns (plus z as its next 
% column if the tracks are 3D). Requires read_gdf.m.

% Written by Nick Ouellette and Doug Kelley. 
% Improved input parsing 28 October 2009.
% 24 February 2010: Now plots are always connected. Also changed "bigfile"
% to "readwhole". 
% Works with 2D or 3D tracks as of 18 October 2011. 

indcol = 1;
xcol = 2;
ycol = 3;
zcol = 4;
%magicnum = 82991; % identifies .gdf files
noisydefault = 1;
readwholedefault = true;
minlengthdefault = 20;

if nargin < 1
    error(['Usage: [nlong,ntracks,lmean,lrms,lmax] = ' mfilename ...
        '(filename,[minlength],[noisy],[readwhole])']);
end
if ~exist('minlength','var') || isempty(minlength)
    minlength = minlengthdefault;
end
if ~exist('noisy','var') || isempty(noisy)
        noisy = noisydefault;
end
if ~exist('readwhole','var') || isempty(readwhole)
    readwhole=readwholedefault;
end

is3D=false;
if readwhole % For smaller .gdf files, read & plot all at once.
    alltracks = filename; %load(filename);%transpose(read_gdf(filename));
    alltracks = alltracks';
    if size(alltracks,1)==6 %should this be 4 (tr, x, y, and z)?
%         disp('Found 6 columns; assuming tracks are 3D.')
        is3D=true;
    end
    [~,ends,~]=unique(alltracks(indcol,:),'legacy');
    begins=circshift(ends,[0,1])+1;
    begins(1)=1;
    lengths=ends-begins+1;
    plotme=lengths>=minlength;
    begins=begins(plotme);
    ends=ends(plotme);

    if noisy
        figure;
        h=axes;
        hold(h, 'on');
        for ii=1:numel(begins)
            ind=begins(ii):ends(ii);
            if is3D
                plot3(alltracks(xcol,ind),alltracks(ycol,ind), ...
                    alltracks(zcol,ind))
            else
                plot(alltracks(xcol,ind),alltracks(ycol,ind))
            end
        end
    end
    
%else % For .gdf files too big for memory, read & plot track-by-track.
%    if noisy
%        figure;
%        axes;
%        hold on;
%    end
%
%    fidin=fopen(filename);
%    lengths=[];
%    header=fread(fidin,6,'int32');
%     if header(1)~=magicnum
%         error(['Sorry, ' filename ' does not appear to be a .gdf file.']);
%     end
%    ncol=header(3); % number of columns of data
%    if ncol==6
%         disp('Found 6 columns; assuming tracks are 3D.')
%        is3D=true;
%    end
%
%    buffer=fread(fidin,ncol,'float32')';
%    while ~feof(fidin) % loop over track index
%        ind=buffer(indcol);
%        fseek(fidin,-4*ncol,'cof'); % back up one sample
%        track=[];
%        while buffer(indcol)==ind % read one whole track
%            buffer=fread(fidin,ncol,'float32')';
%            track=[track;buffer];
%            if feof(fidin)
%                break
%            end
%        end
%        if ~feof(fidin)
%            track(end,:)=[]; % remove last sample, which comes from next track
%        end
%        lengths=[lengths;length(track)];
%        if noisy
%            if length(track)>=minlength
%                if is3D
%                    plot3(track(:,xcol),track(:,ycol),track(:,zcol))
%                else
%                    plot(track(:,xcol),track(:,ycol));
%                end
%            end
%        end
%    end % while ~feof(fidin)
%end % if readwhole
 
nlong=sum(lengths>=minlength); %number of tracks longer than minlength
ntracks=numel(lengths); %total number of tracks
lmean=mean(lengths); %mean track length 
lrms=sqrt(mean(lengths.^2)); %root-mean-square track length
lmax=max(lengths); %maximum track length

if noisy 
    axis tight
    set(gca,'dataaspectratio',[1 1 1])
    fig_title = "Tracks longer than "+int2str(minlength)+" frames";
    title(fig_title,'interpreter','none')
    if is3D
        view(3)
    end
end

end % function longtracks

