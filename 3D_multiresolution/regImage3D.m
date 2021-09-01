function [ f,count,I0, Vx_new, Vy_new, Vz_new] = regImage3D(I0_ori,filename_prefix,TOL,level,nlevel,f,count, I00,I0,I1, Vxl, Vyl, Vzl,gamma,max_Iteration,PlotAll,UseGaussian,SaveVTK,compute_DiceSimilarity,fid_out)
% regImage3D register two given image (I0 and I1) by defoemting the image
% in the normal direction and update the image by B-spline composition

% input: I0_ori = Original moving image
%        filename_prefix: filename_prefix for saving the vtk file
%        TOL = stopping tolerance for the registration process
%        level = multiresolution level
%        nlevel = number of multiresolution level
%        f = frame plotting slices of the images during registration process
%        count =  accumulated number of frame
%        I00 = original moving image scaled according to multiresolution level
%        I0 = moving image scaled acoording to multiresolution level
%        I1 = fixed image scaled acoording to multiresolution level
%        Vxl, Vyl,Vzl = transformation field scaled according to multiresolution level
%        gamma = regularize parameter used to control diffeomorphic transformation
%        max_Iteration = maxi number of iteration
%        PlotAll = option used to set whether to plot the image slices or not
%        UseGaussian = option used to set whether to use gaussian filter to smooth the tranformation field
%        SaveVTK =option used to set whether to save moving image as vtk file or not
%        compute_DiceSimilarity = option used to set whether to compute dice similarity or not
%        fid_out = file location used to save dice similarity 

% output : f = saved frame
%         count = accumulated number of frame
%         I0 = registed moving image
%         Vx_new, Vy_new, Vz_new = resulting transformation field.

[ Bsplinekernel ] = BsplineKernel3D;
[ dxKernel, dyKernel, dzKernel ] = DiffBsplineKernel3D;

ny = size(I0,1);
nx = size(I0,2);
nz = size(I0,3);

I0x=0*meshgrid(1:nz,1:ny);
I0y=0*meshgrid(1:nz,1:nx);
I1x=0*meshgrid(1:nz,1:ny);
I1y=0*meshgrid(1:nz,1:nx);
[Vx_ident, Vy_ident, Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));

ny_ori = size( I0_ori,1);
nx_ori = size( I0_ori,2);
nz_ori = size( I0_ori,3);
[plot_x, plot_y, plot_z] = meshgrid((1:size(I0_ori,2)),(1:size(I0_ori,1)),(1:size(I0_ori,3)));

position_slice_x=round(nx/2);
position_slice_y=round(ny/2);
position_slice_z=round(nz/2);

if compute_DiceSimilarity && PlotAll
    for i=1:nz
        I0x(:,i)=I00(:,position_slice_x,i);
    end
    
    I0xx=I0x';
    I0xx=flipud(I0xx);
    subplot(3,4,1)
    imagesc(I0xx)
    title(['Slice of I0 at x= ', num2str(position_slice_x)])
    axis image
    colormap('gray')
    colorbar
    drawnow
    
    %slice at y position
    for i=1:nz
        I0y(:,i)=I00(position_slice_y,:,i);
    end
    
    I0yy=I0y';
    I0yy=flipud(I0yy);
    subplot(3,4,5)
    imagesc(I0yy)
    title(['Slice of I0 at y= ', num2str(position_slice_y)])
    axis image
    colormap('gray')
    colorbar
    drawnow
    
    %slice at z position
    subplot(3,4,9)
    imagesc(flipud(I00(:,:,position_slice_z)'))
    title(['Slice of I0 at z= ', num2str(position_slice_z)])
    axis image
    colormap('gray')
    colorbar
    drawnow
    
    for i=1:nz
        I1x(:,i)=I1(:,position_slice_x,i);
    end
    
    I1x=I1x';
    I1x=flipud(I1x);
    subplot(3,4,3)
    imagesc(I1x)
    title(['Slice of I1 at x= ', num2str(position_slice_x)])
    axis image
    colormap('gray')
    colorbar
    drawnow
    
    %slice at y position
    for i=1:nz
        I1y(:,i)=I00(position_slice_y,:,i);
    end
    
    I1yy=I1y';
    I1yy=flipud(I1yy);
    subplot(3,4,7)
    imagesc(I1yy)
    title(['Slice of I1 at y= ', num2str(position_slice_y)])
    axis image
    colormap('gray')
    colorbar
    drawnow
    
    %slice at z position
    subplot(3,4,11)
    imagesc(flipud(I1(:,:,position_slice_z)'))
    title(['Slice of I1 at z= ', num2str(position_slice_z)])
    axis image
    colormap('gray')
    colorbar
    drawnow
    
end


I00 = I00(:)';
coef_I0= img2coef3D(I00,nx,ny,nz);
coef= img2coef3D(I0(:)',nx,ny,nz);

%small number used to avoid division by zero in evaluating normal
%direction.
eps = 0.001;

Vxl = Vxl(:)';
Vyl = Vyl(:)';
Vzl = Vzl(:)';

coef_x = img2coef3D(Vxl,nx,ny,nz);
coef_y = img2coef3D(Vyl,nx,ny,nz);
coef_z = img2coef3D(Vzl,nx,ny,nz);


residualHistory=ones(1,6);

TOL_Average=1*10^-(5-level);

residual_initial = norm3d(I1-I0);

residual=norm3d(I1-I0);
relative_residual=residual/residual_initial;

SR(1)=relative_residual;
SRcount=2;
[ I1label ] = ConvertIntensity( I1 );
LabelSize=12;
for k=1:max_Iteration
    
    if relative_residual>TOL
        
        %compute partial derivative
        dIdx_temp = imfilter(coef, dxKernel);
        dIdy_temp = imfilter(coef, dyKernel);
        dIdz_temp = imfilter(coef, dzKernel);
        
        dIdx = dIdx_temp(2:end-2, 2:end-2,2:end-2);
        dIdy = dIdy_temp(2:end-2, 2:end-2,2:end-2);
        dIdz = dIdz_temp(2:end-2, 2:end-2,2:end-2);
        
        N=sqrt((dIdx).^2+(dIdy).^2+(dIdz).^2+eps);
        
        time_step = 1;
        
        %compute transformation field
        Vx_temp = time_step*(I1-I0).*(dIdx)./(N.^2+gamma*(I1-I0).^2);
        Vy_temp = time_step*(I1-I0).*(dIdy)./(N.^2+gamma*(I1-I0).^2);
        Vz_temp = time_step*(I1-I0).*(dIdz)./(N.^2+gamma*(I1-I0).^2);
        
        if UseGaussian
            Vx_temp= smooth3(Vx_temp,'gaussian',[5 5 5]);
            Vy_temp= smooth3(Vy_temp,'gaussian',[5 5 5]);
            Vz_temp= smooth3(Vz_temp,'gaussian',[5 5 5]);
        end
        
        VX=Vx_temp+Vx_ident;
        VY= -Vy_temp+Vy_ident;
        VZ=Vz_temp+Vz_ident;
        
        VX=VX(:)';
        VY=VY(:)';
        VZ=VZ(:)';
        
        CX=coef_x(:)';
        CY=coef_y(:)';
        CZ=coef_z(:)';
        
        %compose transformation field
        [ Vx_new, Vy_new, Vz_new] = BsplineCompose3D( nx,ny,nz,VX, VY, VZ, CX, CY, CZ);
        
        % save transfoamtion it interms of coefficient
        temp_coef_x = imfilter(Vx_new, Bsplinekernel);
        temp_coef_y = imfilter(Vy_new, Bsplinekernel);
        temp_coef_z = imfilter(Vz_new, Bsplinekernel);
        
        coef_x(4:ny,4:nx,4:nz) = temp_coef_x(2:end-2,2:end-2,2:end-2);
        coef_y(4:ny,4:nx,4:nz) = temp_coef_y(2:end-2,2:end-2,2:end-2);
        coef_z(4:ny,4:nx,4:nz) = temp_coef_z(2:end-2,2:end-2,2:end-2);
        
        VX_new=Vx_new(:)';
        VY_new=Vy_new(:)';
        VZ_new=Vz_new(:)';
        
        CI0=coef_I0(:)';
        
        %compose transformation field with the image
        [ I0 ] = BsplineComposeImage3D(nx,ny,nz, VX_new, VY_new, VZ_new, CI0);
        
        if compute_DiceSimilarity

            [ I0label ] = ConvertIntensity( I0 );
            
            II0=I0label(:)';
            II1=I1label(:)';
            
            [ Dice ] = DiceSimilarity( nx, ny, nz, II0, II1 ,LabelSize);
            
            fprintf(fid_out, '%d, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %d\r\n',nlevel-level+1, k, Dice(1), Dice(2), Dice(3), Dice(4), Dice(5), Dice(6),Dice(7), Dice(8), Dice(9), Dice(10), Dice(11), Dice(12));
          
            
        end

        difference_after_deformation=I1-I0;
        residual = norm3d(difference_after_deformation);
        relative_residual = residual/residual_initial;
        disp(['relative_residual = ',num2str(relative_residual)])
        
        temp_coef=imfilter(I0, Bsplinekernel );
        coef(4:ny,4:nx,4:nz)=temp_coef(2:end-2,2:end-2,2:end-2);
        
        if level==1
            I0_vtk=I0;
        else
            I0_vtk = resize3DmatrixScaleIntensity(nx_ori,ny_ori,nz_ori,I0);
        end
        
        if SaveVTK
            filename = [filename_prefix, num2str(count),'.vtk'];
            vtkwrite(filename,'structured_grid',plot_x,plot_y,plot_z,'scalars','Intensity',I0_vtk);
        end
        
        if PlotAll && compute_DiceSimilarity==0
            SR(SRcount)=relative_residual;
            SR_plot=SR(1:SRcount);
            T=linspace(1,SRcount,SRcount);
            th=subplot(2,2,2);
            plot(T,SR_plot)
            hold on
            title('Graph of relative residual')
            xlabel('Iteration')
            ylabel('Similarity ratio')
            SRcount=SRcount+1;
            
            %slice at x position
            for i=1:nz
                I0x(:,i)=I0(:,position_slice_x,i);
            end
            I0xx=I0x';
            I0xx=flipud(I0xx);
            subplot(2,2,1)
            imagesc(I0xx)
            title({['moving image at step:', num2str(k)];['at x= ', num2str(position_slice_x)]})
            axis image
            colormap('gray')
            colorbar
            drawnow
            
            for i=1:nz
                I0y(:,i)=I0(position_slice_y,:,i);
            end
            
            I0yy=I0y';
            I0yy=flipud(I0yy);
            subplot(2,2,3)
            imagesc(I0yy)
            title({['moving image at y=', num2str(position_slice_y)];['at level= ', num2str(nlevel-level+1)]})
            axis image
            colormap('gray')
            colorbar
            drawnow
            
            subplot(2,2,4)
            imagesc(flipud(I0(:,:,position_slice_z)'))
            title(['moving Image at z= ', num2str(position_slice_z)])
            axis image
            colormap('gray')
            colorbar
            drawnow
            f(count) = getframe(gcf, [0 0 1200 800]);
        end
        
        if compute_DiceSimilarity && PlotAll
            for i=1:nz
                I0x(:,i)=I0(:,position_slice_x,i);
            end
            
            I0xx=I0x';
            I0xx=flipud(I0xx);
            subplot(3,4,2)
            imagesc(I0xx)
            title({['Image registration at step:', num2str(k)];['Slice of deformed I0 at x= ', num2str(position_slice_x)]})
            axis image
            colormap('gray')
            colorbar
            drawnow
            
            %slice at y position
            for i=1:nz
                I0y(:,i)=I0(position_slice_y,:,i);
            end
            
            I0yy=I0y';
            I0yy=flipud(I0yy);
            subplot(3,4,6)
            imagesc(I0yy)
            title({['Slice of deformed I0 at y= ', num2str(position_slice_y)];['at level= ', num2str(nlevel-level+1)]})
            axis image
            colormap('gray')
            colorbar
            drawnow
            
            %slice at z position
            subplot(3,4,10)
            imagesc(flipud(I0(:,:,position_slice_z)'))
            title(['Slice of deformed I0 at z= ', num2str(position_slice_z)])
            axis image
            colormap('gray')
            colorbar
            drawnow
            
            SR(SRcount)=relative_residual;
            SR_plot=SR(1:SRcount);
            T=linspace(1,SRcount,SRcount);
            th=subplot(3,4,4);
            plot(T,SR_plot)
            hold on
            title('Graph of relative residual')
            xlabel('Iteration')
            ylabel('Similarity ratio')
            SRcount=SRcount+1;
            
            f(count) = getframe(gcf, [0 0 1200 800]);
            
        end
        
        residualHistory(1:5) = residualHistory(2:6);
        residualHistory(6) = relative_residual;
        avgresidualHistory1=(residualHistory(1)+residualHistory(2)+residualHistory(3))/3;
        avgresidualHistory2=(residualHistory(4)+residualHistory(5)+residualHistory(6))/3;
        
        if avgresidualHistory1 -avgresidualHistory2 < TOL_Average
            break
        end
        
        count=count+1;
        
    else
        
        break;
    end
    toc
end

end

