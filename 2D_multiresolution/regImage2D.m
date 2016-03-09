function [f,I0,count,Vx, Vy, Vx_inv, Vy_inv,FgridX, FgridY,In_gridX, In_gridY] = regImage2D( scale,TOL_Average,level,nlevel, f,count,I0_ori, I0,I1, Vxl, Vyl, Vxl_inv, Vyl_inv,M_filt,Fx, Fy,max_Iteration,PlotAll,UseGaussian,TOL,gamma)
%Input:  scale               = scale used to downsample image for registration at coarser resolution
%        TOL_Average         = stopping tolerance for average of relative residual
%        level               = multiresolution level
%        nlevel              = number of multiresolution level
%        f                   = frame for plotting images
%        count               = number  of frames
%        I0_ori              = moving image with default size
%        I0                  = scaled moving image
%        I1                  = scaled fixed image
%        Vxl, Vy1            = Identity transformation in x direction and y direction
%        Vxl_inv, Vyl_inv    = Identity inverse transformation in x direction and y direction
%        M_filt              = Bspline kernel
%        Fx & Fy             = discrete kernel for computing partial derivative with respect to x and y direction respectively.  
%        max_Iteration       = maximum number of iteration
%        PlotAll             = option whether to plot images throughout the registration process or not
%        UseGaussian         = option whether to use Gaussion filter to smooth the transformation field 
%        TOL                 = stopping tolerance or not
%        gamma               = regularize parameter for diffeomorphic registration. 

% Output :  f                       = saved frame
%           I0                      = registered moving image
%           Vx & Vy                 = resulting transformation field in x and y direction respectively 
%           Vx_inv & Vy_inv         = resulting inverse transformation field in x and y direction respectively
%           FgridX & FgridY         = resulting forward defomation grid in x and y direction respectively
%           In_gridX & In_gridY     = resulting inverse deformation grid in x and y direction respectively

% degree of B-spline
p=3;
q=3;

if PlotAll
    subplot(2,4,1)
    imagesc(I0)
    title('Source image')
    colormap('gray')
    %colorbar
    axis image
    drawnow
    
    subplot(2,4,2)
    imagesc(I1)
    title('Target image')
    axis image
    colormap('gray')
    drawnow
end

%gaussian filter
sigma_par = round(10*scale);
smooth_par = 6*sigma_par;
Hsmooth=fspecial('gaussian',[smooth_par, smooth_par], sigma_par);

ny = size(I0,1);
nx = size(I0,2);

II0=I0(:)';
II0_ori=I0_ori(:)';

coef= img2coef2D(nx,ny,II0);
coef_I0= img2coef2D(nx,ny,II0_ori);
CI0 = coef_I0(:)';


%compute differecnt between moving image and fixed image
residual_initial = norm(I1-I0);

residual=norm(I1-I0);
relative_residual=residual/residual_initial;
I_diff=I1-I0;
if PlotAll
    subplot(2,4,3)
    imagesc(I_diff)
    title(['Current Diff. Rel. resid=',num2str(relative_residual)])
    axis image
    colormap('gray')
    drawnow
end
alpha = 0.001;

%initialize the deformed grid
ngx = 40;
ngy = 40;
VGx = zeros(ngy,ngx);
VGy = zeros(ngy,ngx);

[FgridX, FgridY, plot_pix_index_column, plot_pix_index_row] = makeGrid(ngx,ngy,nx,ny);
gridX0 = FgridX;
gridY0 = FgridY;

if PlotAll
    th3=subplot(2,4,6);
    cla(th3)
    plotGrid(FgridX, FgridY)
    axis image
    title('Fwd. deformed grid')
    drawnow
    hold off
end

VXL=Vxl(:)';
VYL=Vyl(:)';
coef_x = img2coef2D(nx,ny, VXL);
coef_y = img2coef2D(nx,ny, VYL);

VXL_inv=Vxl_inv(:)';
VYL_inv=Vyl_inv(:)';

coef_x_inv = img2coef2D(nx, ny, VXL_inv);
coef_y_inv = img2coef2D(nx, ny, VYL_inv);

residualHistory=10*ones(1,6);


[Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5), (0.5:ny-0.5));

SR(1)=relative_residual;
SRcount=2;


for k=1:max_Iteration
    
    if relative_residual>TOL
        disp(['Iteration = ' num2str(k)]) ;
        dIdx_temp = imfilter(coef, Fx);
        dIdy_temp = imfilter(coef, Fy);
        
        dIdx = dIdx_temp(2:end-2, 2:end-2);
        dIdy = dIdy_temp(2:end-2, 2:end-2);
        
        N=sqrt((dIdx).^2+(dIdy).^2+alpha);
        
        time_step = 1;

        Vx_temp = time_step*(I1-I0).*(dIdx)./(N.^2+gamma*(I1-I0).^2);
        Vy_temp = time_step*(I1-I0).*(dIdy)./(N.^2+gamma*(I1-I0).^2);
        
        if UseGaussian
            Vx_temp=imfilter(Vx_temp, Hsmooth);
            Vy_temp=imfilter(Vy_temp, Hsmooth);
        end
        
        Vx_old = Vx_temp+Vx_ident;
        Vy_old = -Vy_temp+Vy_ident;
        
        Vx_old = reshape(Vx_old, 1, nx*ny);
        Vy_old = reshape(Vy_old, 1, nx*ny);
        
        coef_x_fun = reshape(coef_x, 1, (nx+p)*(ny+q));
        coef_y_fun = reshape(coef_y, 1, (nx+p)*(ny+q));
        
        Vx_old_inv = -Vx_temp+Vx_ident;
        Vy_old_inv = Vy_temp+Vy_ident;
        
        Vx_old_inv = reshape(Vx_old_inv, 1, nx*ny);
        Vy_old_inv = reshape(Vy_old_inv, 1, nx*ny);
        
        coef_x_inv_fun = reshape(coef_x_inv, 1, (nx+p)*(ny+q));
        coef_y_inv_fun = reshape(coef_y_inv, 1, (nx+p)*(ny+q));
        
        [ Vx_new, Vy_new] = BsplineCompose2D( Vx_old, Vy_old, coef_x_fun, coef_y_fun, nx, ny);
        [ Vx_new_inv, Vy_new_inv] = BsplineCompose2D( Vx_old_inv, Vy_old_inv, coef_x_inv_fun, coef_y_inv_fun, nx, ny);
        
        temp_coef_x = imfilter(Vx_new, M_filt);
        temp_coef_y = imfilter(Vy_new, M_filt);
        
        coef_x(4:ny,4:nx) = temp_coef_x(2:end-2,2:end-2);
        coef_y(4:ny,4:nx) = temp_coef_y(2:end-2,2:end-2);
        
        temp_coef_x_inv = imfilter(Vx_new_inv, M_filt);
        temp_coef_y_inv = imfilter(Vy_new_inv, M_filt);
        
        coef_x_inv(4:ny,4:nx) = temp_coef_x_inv(2:end-2,2:end-2);
        coef_y_inv(4:ny,4:nx) = temp_coef_y_inv(2:end-2,2:end-2);
        
        Vx = Vx_new - Vx_ident;
        Vy = Vy_new - Vy_ident;
        
        Vx_inv = Vx_new_inv - Vx_ident;
        Vy_inv = Vy_new_inv - Vy_ident;
        
        Vy=-Vy;
        Vy_inv = -Vy_inv;
        
        Vx_new = reshape(Vx_new, 1, nx*ny);
        Vy_new = reshape(Vy_new, 1, nx*ny);
        
        I0 = BsplineComposeImage2D(Vx_new, Vy_new, CI0, nx, ny);
        
        In_gridX=zeros(size(FgridX));
        In_gridY=zeros(size(FgridY));
        
        for ii=1:ngx
            for jj=1:ngy
                
                j = gridY0(jj,ii);
                i = gridX0(jj,ii);
                
                VGx(jj,ii) = Vx(j,i);
                VGy(jj,ii) = Vy(j,i);
                
                In_gridX(jj,ii) = round(i+VGx(jj,ii));
                In_gridY(jj,ii) = round(j-VGy(jj,ii));
                
            end
        end
        
        if PlotAll
            th2=subplot(2,4,4);
            cla(th2);
            plotGrid(In_gridX, In_gridY)
            axis image
            title(['Inverse deformed grid at step: ',num2str(k) ])
            colormap('gray')
            set(gca,'YDir','reverse');
            drawnow
            hold off
        end
        
        for ii=1:ngx
            for jj=1:ngy
                
                j = gridY0(jj,ii);
                i = gridX0(jj,ii);
                
                VGx(jj,ii) = Vx_inv(j,i);
                VGy(jj,ii) = Vy_inv(j,i);
                
                FgridX(jj,ii) = round(i+VGx(jj,ii));
                FgridY(jj,ii) = round(j-VGy(jj,ii));
                
            end
        end
        
        if PlotAll
            subplot(2,4,5);
            imagesc(I0)
            title(['Current image at level: ',num2str(nlevel-level+1) ])
            axis image
            colormap('gray')
            drawnow
            
            th1=subplot(2,4,6);
            cla(th1);
            plotGrid(FgridX, FgridY)
            axis image
            title(['Deformed grid at step: ',num2str(k) ])
            colormap('gray')
            set(gca,'YDir','reverse');
            drawnow
            hold off
            
            quiver_matrix_x = Vx(plot_pix_index_row, plot_pix_index_column);
            quiver_matrix_y = Vy(plot_pix_index_row, plot_pix_index_column);
            th5=subplot(2,4,7);
            cla(th5);
            quiver(-quiver_matrix_x(end:-1:1,:),-quiver_matrix_y(end:-1:1,:),2)
            title(['Transformation field at step ', num2str(k)])
            axis tight
            drawnow
        end
        
        difference_after_deformation=I1-I0 ;
        residual = norm(difference_after_deformation);
        relative_residual = residual/residual_initial;
        disp(['Relative_residual = ' num2str(relative_residual)]) ;
        
        if PlotAll
            %plot relative residual
            SR(SRcount)=relative_residual;
            SR_plot=SR(1:SRcount);
            T=linspace(1,SRcount,SRcount);
            th=subplot(2,4,8);
            plot(T,SR_plot)
            hold on
            title('Graph of relative residual')
            xlabel('Iteration')
            ylabel('Similarity ratio')
            SRcount=SRcount+1;

            th4=subplot(2,4,3);
            cla(th4);
            imagesc(I1-I0)
            title(['Current Diff. Rel. resid=',num2str(relative_residual)])
            axis image
            colormap('gray')
            drawnow
        end
        
        temp_coef=imfilter(I0, M_filt);
        coef(4:ny,4:nx)=temp_coef(2:end-2,2:end-2);
        
        if PlotAll
            f(count) = getframe(gcf, [0 0 1200 800]);
            count=count+1;
        end
        
        residualHistory(1:5) = residualHistory(2:6);
        residualHistory(6) = relative_residual;
        avgresidualHistory1=(residualHistory(1)+residualHistory(2)+residualHistory(3))/3;
        avgresidualHistory2=(residualHistory(4)+residualHistory(5)+residualHistory(6))/3;
        
        if avgresidualHistory1 -avgresidualHistory2 < TOL_Average
            break
        end
        
    else
        
        break;
    end
    
    toc
    
end
end

