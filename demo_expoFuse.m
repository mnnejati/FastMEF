%----------------------------------------------------
%  A demo of fast multi-exposure image fusion method
%
%----------------------------------------------------
clc,clear all,close all;
addpath('utils');

%========< load input image sequence >==========
imgName = 'Tower';  % Tower , Kluki, ...
imgPath = 'images\';
files = dir(fullfile(imgPath,imgName,'\','*.png'));
N = length(files);
for i = 1:N
    Ii = imread(['images\' imgName '\' files(i).name]);
    Io(:,:,:,i) = im2double(Ii);
end


%========< multi-exposure image fusion >==========
wRad = 12; % radius of the filter

tic;
FI = fastExpoFuse(Io, wRad);
toc;
figure,imshow(FI)
