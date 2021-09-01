function [I0,Vx,Vy] = MultiresolutionRegistration2D_silent(I0,I1,nlevel,max_Iteration,PlotAll,SaveVideo,filename,UseGaussian,TOL,gamma)

[~,I0,Vx,Vy] = evalc(['MultiresolutionRegistration2D(I0,I1,',...
    'nlevel,max_Iteration,PlotAll,SaveVideo,filename,UseGaussian,TOL,gamma)']);

end