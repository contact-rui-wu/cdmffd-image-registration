%2D multiresolution image registration usign level set and B-spline compose
%Lena example
close all
clear all
tic
addpaths

filename_i0 = 'lena2.png'; I0=double(imread(filename_i0));
filename_i1 = 'lena1.png'; I1=double(imread(filename_i1));
filename='Lena.avi';

nlevel=3;
max_Iteration=1000;
TOL=0.05;
gamma=1.54;
PlotAll = 1;
SaveVideo = 0;
UseGaussian=0;

[I0] = MultiresolutionRegistration2D( I0,I1,nlevel,max_Iteration,PlotAll,SaveVideo,filename,UseGaussian,TOL,gamma);
