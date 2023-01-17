function [I0,Vx,Vy,Vz,output] = MultiresolutionRegistration3D_silent(I0,I1,nlevel,max_Iteration,PlotAll,SaveVideo,SaveVTK,compute_DiceSimilarity,TOL,UseGaussian,gamma)

[~,I0,Vx,Vy,Vz,output] = evalc(['MultiresolutionRegistration3D(I0,I1,'...
    'nlevel,max_Iteration,PlotAll,SaveVideo,SaveVTK,compute_DiceSimilarity,'...
    'TOL,UseGaussian,gamma)']);

end