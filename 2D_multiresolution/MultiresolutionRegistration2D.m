function [I0,Vx,Vy,output] = MultiresolutionRegistration2D( I0,I1,nlevel,max_Iteration,PlotAll,SaveVideo,filename,UseGaussian,TOL,gamma)
%regImage_multiresolution2D register two given images from coarse to fine
%level using level set and B-spline composition
%Input:  I0               = moving image
%        I1               = fixed image
%        nlevel           = multiresolution level
%        max_Iteration    = maximum Iteration in each level
%        PlotAll          = option used to set whether to plot images during the registration process or not.
%        SaveVideo        = option used to set whether to save the plotted frame as video or not. 
%        filename         = filename of the video
%        UseGaussian      = option used to set whether to use a Gaussian filter to smooth the transformation field or not. 
%        TOL              = stopping tolerance
%        gamma            = regularize parameter for diffeomorphic registration. 

%output: I0 = registration result


I00=I0;
I11=I1;

%Basis value computed by substituting corresponding midpoint
M = [1/24, 11/24, 11/24, 1/24];

%Bspline kernel
M_filt = M'*M;

% discrete kernel for computing partial derivative with respect
% to x
Fx =[1/48, 23/48,  23/48, 1/48]'*[-1/8 -5/8, 5/8 1/8];
% discrete kernel for computing partial derivative with respect
% to y
Fy =[1/8 5/8,  -5/8 -1/8]'*[1/48, 23/48,  23/48, 1/48];

count=1;
if PlotAll
figure('position',[100 100 1200 800]);
f(1) = getframe(gcf, [0 0 1200 800]);
else
f=[];    
end
levelcount=1;

tic
for level=nlevel:-1:1
    disp(['Register level: ' num2str(nlevel-level+1) '...']);
      
    disp('Scaling image to the corresponding size of current multiresolution level...');
    % downsample image
    scale = 2^-(level-1);
    I1 = imresize(I11,scale);
    I0 = imresize(I0,scale);
    I0_ori=imresize(I00,scale);
    
    if level==nlevel
        
        nx=size(I0,2);
        ny=size(I0,1);
        [Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5), (0.5:ny-0.5));
        Vxl = Vx_ident;
        Vyl = Vy_ident;
        
        Vxl_inv = Vx_ident;
        Vyl_inv = Vy_ident;
        
    else
        
        nx=size(I0,2);
        ny=size(I0,1);
        [Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5), (0.5:ny-0.5));
        
        disp('Scaling the transformation field to the corresponding size of current multiresolution level...');
        Vxl = imresize(Vx*scale,scale);
        Vyl = imresize(Vy*scale,scale);
        
        Vxl_inv = imresize(Vx_inv*scale,scale);
        Vyl_inv = imresize(Vy_inv*scale,scale);
        
        Vxl=Vxl+Vx_ident;
        Vyl=-Vyl+Vy_ident;
        
        Vxl_inv=Vxl_inv+Vx_ident;
        Vyl_inv=-Vyl_inv+Vy_ident;
        
        II0_ori=I0_ori(:)';
        
        disp('Computing coefficient of level set representation of the scaled image...');
        coef_I0= img2coef2D(nx,ny, II0_ori);
        
        CI0 = coef_I0(:)';
        
        VXL=Vxl(:)';
        VYL=Vyl(:)';
        disp('Composing transformation field from previous level to the image...');
        I0 = BsplineComposeImage2D(VXL, VYL, CI0, nx, ny);
    end
    
    TOL_Average=10^-((4+nlevel)-level);
    % registration
    [f,I0,count, Vx, Vy, Vx_inv, Vy_inv,FgridX, FgridY,In_gridX, In_gridY] = regImage2D(scale,TOL_Average,level,nlevel, f,count,I0_ori, I0,I1, Vxl, Vyl, Vxl_inv, Vyl_inv,M_filt,Fx, Fy,max_Iteration,PlotAll,UseGaussian,TOL,gamma);
    
    output(level).Vx = Vx;
    output(level).Vy = Vy;
    output(level).FgridX = FgridX;
    output(level).FgridY = FgridY;
    output(level).regIm = I0;
    output(level).scale = scale;
    
    disp('Scaling the transformation field back to its original size...');
    % upsample
    Vx = imresize(Vx/scale,size(I00));
    Vy = imresize(Vy/scale,size(I00));
    
    Vx_inv = imresize(Vx_inv/scale,size(I00));
    Vy_inv = imresize(Vy_inv/scale,size(I00));
    
    disp('Scaling the image back to its original size...');
    I0=imresize(I0/scale,size(I00));
    
    levelcount=levelcount+1;
    
end

InitialResidual=norm(I11-I00);
FinalResidual=norm(I11-I0);
Final_relative_residual=FinalResidual/InitialResidual;
disp(['Final Relative Residual =',num2str(Final_relative_residual)])


figure('position',[100 100 1200 800]);
subplot(2,3,1)
imagesc(I0_ori)
title('Source image')
colormap('gray')
axis image
drawnow

subplot(2,3,2)
imagesc(I1)
title('Target image')
axis image
colormap('gray')
drawnow

subplot(2,3,3);
imagesc(I0)
title('Registered Moving image')
axis image
colormap('gray')
drawnow

subplot(2,3,4)
plotGrid(FgridX, FgridY)
axis image
title('Fwd. deformed grid')
drawnow
hold off

subplot(2,3,5)
plotGrid(In_gridX, In_gridY)
axis image
title('Inverse deformed grid')
colormap('gray')
set(gca,'YDir','reverse');
drawnow
hold off

subplot(2,3,6);
imagesc(I1-I0)
title(['Diff. & Final Rel. resid=',num2str(Final_relative_residual)])
axis image
colormap('gray')
drawnow


if SaveVideo
frameRate = 15;
vidObj = VideoWriter(filename);
vidObj.FrameRate = frameRate;
open(vidObj)
writeVideo(vidObj, f)
close(vidObj)
end


end

