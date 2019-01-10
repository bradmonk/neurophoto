function [PIX] = getIMGpaths()


% Open image selection dialogue box
[fname,fpath] = uigetfile({'*.tif*';'*.TIF*'});
imgfullpath = [fpath,fname];

% Select image to import
PIX.path = imgfullpath;


% Use iminfo() to gather information about image file
Info = imfinfo(PIX.path);


% Turn the struct into a table
PIX.info = struct2table(Info);


% Count the number of images in the stack
PIX.count  = height(PIX.info);



end