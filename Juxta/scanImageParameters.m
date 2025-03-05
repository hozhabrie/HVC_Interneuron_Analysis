clear all; close all; clc;

%first download scanimage tiffreader folder

%then import this line and set the path
import ScanImageTiffReader.ScanImageTiffReader;
path = 'C:\Users\User\Dropbox (NYU Langone Health)\Int Juxta\032819_dlx32\';

% list the names of the tiffs of the RAW data in the folder
%cell = {'eh9c1'; 'eh16c1'; 'eh26c3'; 'eh33c2'; 'eh36c4'; 'eh51c4'};
cell = {'cell1_zoom4_2002'}
%for every cell...
for i = 1:length(cell)
    %read the tiff as an object
    obj = ScanImageTiffReader([path cell{i} '.tif']);
    
    %obtain the metadata embedded in the tiff 
    meta = obj.descriptions{1};
    
    %set the keywords you want to find
    param = 'state.acq.zoomFactor';
    
    %find where the parameter is in the long string
    loc = strfind(meta,param);
    
    %use this to check your work
    %meta(loc:loc+length(param)+1);
    
    % get the zoom factor (or whatever the param is)
    zoom = meta(loc+length(param)+1);
    
    % get the distance in fov: 608 um x 608 um is zoom 1
    fov = 608/str2num(zoom);
    
    %display results
    sprintf([string(cell{i}) + '\nZoom = ' + zoom + '\nFOV = ' + string(fov) + ' microns'])
end

% then open AVERAGED tiff on Fiji
% Analyze > set scale 
% Dist in pixels = 1024(or whatever the pixels say on the tiff)
% Known distance = fov 
% Unit of Length = microns

% Analyze > Tools > Scalebar
% Width 10  
% Font Size 30
% idc the values just make it obvious and consistent 
% make sure to save the new image with scale bar
