%3D multiresolution image registration usign level set and B-spline compose
%sphere to sun-like shape example
close all
clear all
tic

load Sphere
load sun_like

nlevel=1;
max_Iteration=100;
gamma = 0.1;
TOL=0.05; 

PlotAll = 1;
SaveVideo = 0;
SaveVtk=0;
compute_DiceSimilarity=0;
UseGaussian=0;

I0 = sphere;
I1 = sun_like;


[I0,Vx,Vy,Vz]=MultiresolutionRegistration3D(I0, I1,nlevel,max_Iteration,PlotAll,SaveVideo,SaveVtk,compute_DiceSimilarity,TOL,UseGaussian,gamma);
