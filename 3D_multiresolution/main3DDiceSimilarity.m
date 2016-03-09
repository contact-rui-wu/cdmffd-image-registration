
close all
clear all
tic

%%%%%%%%%%%%%%%%%%%%%%%%%

filenameI0='subject04_crisp_v.rawb';
filenameI1='subject05_crisp_v.rawb';

[BrainI0] = readrawb(filenameI0);
[BrainI1] = readrawb(filenameI1);

scale=1/2;
nx=size(BrainI0,2);
ny=size(BrainI0,1);
nz=size(BrainI0,3);
new_size_x=round(nx*scale);
new_size_y=round(ny*scale);
new_size_z=round(nz*scale);

BrainI0 = resize3Dmatrix(new_size_x,new_size_y,new_size_z,BrainI0);
BrainI1 = resize3Dmatrix(new_size_x,new_size_y,new_size_z,BrainI1);

BrainI0=(BrainI0-min(BrainI0(:)))/(max(BrainI0(:))-min(BrainI0(:)))*11;
BrainI0=round(BrainI0);

BrainI1=(BrainI1-min(BrainI1(:)))/(max(BrainI1(:))-min(BrainI1(:)))*11;
BrainI1=round(BrainI1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert image to 0-255
I0=BrainI1*23.1818;
I1=BrainI0*23.1818;

nlevel=3;
max_Iteration=1000;
gamma = 1.54;
PlotAll = 1;
TOL=0.05;

SaveVideo = 0;
SaveVtk = 0;
compute_DiceSimilarity=1;
UseGaussian=0;

[I0]=MultiresolutionRegistration3D(I0, I1,nlevel,max_Iteration,PlotAll,SaveVideo,SaveVtk,compute_DiceSimilarity,TOL,UseGaussian,gamma);
