%2D multiresolution image registration usign level set and B-spline compose
%C-shape example
close all
clear all
tic
addpaths

filename_i0 = 'whitecircle.png'; I0=double(rgb2gray(imread(filename_i0)));
filename_i1 = 'whitec.png'; I1=double(rgb2gray(imread(filename_i1)));
filename='C_shape.avi';

nlevel=3;
max_Iteration=1000;
TOL=0.05;
gamma=1.54;
PlotAll = 1;
SaveVideo = 0;
UseGaussian=0;

[I0] = MultiresolutionRegistration2D( I0,I1,nlevel,max_Iteration,PlotAll,SaveVideo,filename,UseGaussian,TOL,gamma);
