function bg=BackgroundImage(inputnames,outputname)
% Usage: bg=BackgroundImage(inputnames,[outputname])
% Given a movie, background_image calculates the mean pixel values over 
% time, returning the result in "bg" and saving it as an image with the 
% filename "outputname". The movie must be saved as a series of image 
% files, an image stack in .tif or .gif format, or an uncompressed .avi 
% file; specify the movie in "inputnames" (e.g., '0*.png' or 'stack.tif' 
% or 'movie.avi'). To specify an output format, include an extension in 
% "outputname"; '.tif' is used by default. This file can be downloaded 
% from http://leviathan.eng.yale.edu/software. Requires 
% read_uncompressed_avi.m for use with .avi movies. 

% Written 6 April 2010 by Doug Kelley. 
% Fixed error catching of output names w/o extension 4 May 2010.
% Specified uncompressed output explicitly 25 May 2010. 
% Allowed for 16-bit images 17 September 2010. 
% Allowed for color images and files in another directory 13 April 2011. 
% Added image size check and color type check 4 May 2011. 
% Allowed for grayscale images with bit depth between 8 and 16, 17 May 2011. 
% Made compatible with uncompressed avi movies (using 
% read_uncompressed_avi.m) 20 October 2011. 
% Made compatible with tiff and gif stacks 13 February 2012. 

% -=- Set defaults -=-----------------------------------------------------
outputnamedefault='background.tif';
formatdefault='.tif'; % for output, if not specified
nocompress={'.hdf','.png','.tif','.tiff'}; % save these formats uncompressed

% -=- Parse inputs -=-----------------------------------------------------
if nargin<1
    error(['Usage: bg = ' mfilename '([inputnames],[outputname])'])
end
if ~exist('outputname','var') || isempty(outputname)
    outputname=outputnamedefault;
end

% -=- Decide whether avi, stack, or images; set up -=---------------------
[prefix,junk,ext]=fileparts(inputnames);
list=dir(inputnames);
if strcmpi(ext,'.avi')
    movtype='avi';
    if isempty(which('read_uncompressed_avi.m')) % check for req'd helper function
        error(['Sorry, reading .avi files requires ' ...
            'read_uncompressed_avi.m.'])
    end
    movinfo=aviinfo(list.name);
    Nf=movinfo.NumFrames;
    bitdepth=log2(movinfo.NumColormapEntries);
    ct='grayscale'; % we'll just assume!
    mark=[]; % set up for later use with read_uncompressed_avi
elseif numel(list)==1 && ( strcmpi(ext,'.tif') || ...
        strcmpi(ext,'.tiff') || strcmpi(ext,'.gif') ) % single file, looks like an image stack
    movtype='stack';
    movinfo=imfinfo(list.name);
    if numel(movinfo)==1
        error([inputnames ': Only one frame found.'])
    end
    Nf=numel(movinfo);
    bitdepth=movinfo(1).BitDepth;
    ct=movinfo(1).ColorType;
else
    movtype='images';
    movinfo=imfinfo(fullfile(prefix,list(1).name));
    Nf=numel(list);
    bitdepth=movinfo.BitDepth;
    ct=movinfo.ColorType;
end
if ~Nf
    error([inputnames ': No frames found.'])
end

% -=- Loop over frames, keeping running sum -=----------------------------
for ii=1:Nf
    switch movtype
        case('avi')
            [d,mark] = read_uncompressed_avi(fullfile(prefix,list.name), ...
                ii,mark);
            d = d.cdata;
        case('stack')
            d=imread(fullfile(prefix,list.name),ii);
        case('images')
            d=imread(fullfile(prefix,list(ii).name));
    end % switch movtype
    if ii==1
        bg0=double(d);
    else
        bg0=bg0+double(d);
    end % if ii==1
end

% -=- Handle color format & bit depth -=----------------------------------
if strcmp(ct,'grayscale')
    if bitdepth==8
        bg0=uint8(round(bg0/Nf));
    elseif bitdepth<=16
        bg0=uint16(round(bg0/Nf));
    else
        error(['Sorry, grayscale images must have bit depth between ' ...
            '8 and 16.'])
    end
else
    disp('Converting color images to grayscale.')
    bg0=uint8(round(mean(bg0,3)/Nf));
end

% -=- Save file and return output if requested -=-------------------------
[junk,junk,ext]=fileparts(outputname);
if isempty(ext) % no extension; will append one
    imwrite(bg0,[outputname formatdefault],'compression','none');
elseif any(strcmp(ext,nocompress))
    imwrite(bg0,outputname,'compression','none');
else
    imwrite(bg0,outputname);
end
if nargout>0
    bg=bg0;
end