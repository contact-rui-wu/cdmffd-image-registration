%3D-multiresolution
function [I0,Vx,Vy,Vz] = MultiresolutionRegistration3D(I0, I1, nlevel,max_Iteration,PlotAll,SaveVideo,SaveVTK,compute_DiceSimilarity,TOL,UseGaussian,gamma)
% main_regImage3D register image Io to I1
% Input: I0 = moving image
%        I1 = fixed image
%        nlevel = number of multi resolution level
%        max_iteration = maximum number of iteration in each level
%        PlotAll = option whether to plot slices of image throughout the
%        registration process or not
%        SaveVideo = Option whether to save Plotted frames as a video or
%        not
%       compute_DiceSimilarity = option whether to compute dice similarity
%       or not
%       TOL = stopping tolErance
%       UseGaussian = option whether to use a Gaussian filter to smooth the
%       transformation field or not.
%       gamma = regularize parameter for diffeomophic transfomation. 
% Output: I0 = registered moving image.

disp('initializing...')
ori_size_x=size(I0,2);
ori_size_y=size(I0,1);
ori_size_z=size(I0,3);

I1_ori=I1;
I0_ori=I0;

if compute_DiceSimilarity
    
    output_file='Dice_similarity.txt';
    fid_out = fopen(output_file,'w');
    fprintf(fid_out, 'Level Iteration  Background  CSF  GrayMatter  WhiteMatter  Fat  Muscle Muscle/Skin  Skull  vessels  AroundFat  DuraMatter  BoneMarrow\r\n');

    [ I0label ] = ConvertIntensity( I0 );
    [ I1label ] = ConvertIntensity( I1 );
    
    LabelSize=12;
    
    nx=size(I0,2);
    ny=size(I0,1);
    nz=size(I0,3);
    
    II0=I0label(:)';
    II1=I1label(:)';
    
    [ Dice ] = DiceSimilarity( nx, ny, nz, II0, II1 ,LabelSize);
    
    fprintf(fid_out, '%d, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f, %3.5f,%d\r\n', 1, 0, Dice(1), Dice(2), Dice(3), Dice(4), Dice(5), Dice(6),Dice(7), Dice(8), Dice(9), Dice(10), Dice(11), Dice(12));
  
else
    fid_out=[];
end

if PlotAll 
    figure('Position',[100,100,1300,900])
    f(1) = getframe(gcf, [0 0 1200 800]);
else
    f=[];
end

count=1;


filename_prefix='RegistrationExample3D';

tic
for level=nlevel:-1:1
    
    disp(['Register level: ' num2str(nlevel-level+1) '...']);
    
    % downsample
    scale = 2^-(level-1);
    
    new_size_x=round(ori_size_x*scale);
    new_size_y=round(ori_size_y*scale);
    new_size_z=round(ori_size_z*scale);
    
     disp('Scaling image to corresponding size of current multiresolution level...');
    I1 = resize3Dmatrix(new_size_x,new_size_y,new_size_z,I1_ori);
    I0 = resize3Dmatrix(new_size_x,new_size_y,new_size_z,I0);
    I00 =resize3Dmatrix(new_size_x,new_size_y,new_size_z,I0_ori);
    
    if level==nlevel
        
        ny = size(I0,1);
        nx = size(I0,2);
        nz = size(I0,3);
        
        [Vx_ident, Vy_ident, Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));
        
        Vxl = Vx_ident;
        Vyl = Vy_ident;
        Vzl = Vz_ident;
        
    else
        
        ny = size(I0,1);
        nx = size(I0,2);
        nz = size(I0,3);
        
        [Vx_ident, Vy_ident, Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));
        
        disp('Scaling the transformation field to corresponding size of current multiresolution level...');
        Vxl= resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vx*scale);
        Vyl= resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vy*scale);
        Vzl= resize3Dmatrix(new_size_x,new_size_y,new_size_z,Vz*scale);
        
        Vxl=Vxl+Vx_ident;
        Vyl=-Vyl+Vy_ident;
        Vzl=Vzl+Vz_ident;
        
        disp('Computing coefficient of level set representation of the scaled image...');
        coef_I0= img2coef3D(I00(:)',size(I00,2),size(I00,1),size(I00,3));
        VXL=Vxl(:)';
        VZL=Vzl(:)';        
        VYL=Vyl(:)';

        
        CI0=coef_I0(:)';
        disp('Composing transformation field from previos level to the image...');
        [ I0 ] = BsplineComposeImage3D( nx,ny, nz,VXL, VYL, VZL, CI0);
        
    end
    
    [f,count,I0, Vx_new, Vy_new, Vz_new] = regImage3D(I0_ori,filename_prefix,TOL,level,nlevel,f,count, I00,I0,I1, Vxl, Vyl, Vzl,gamma,max_Iteration,PlotAll,UseGaussian,SaveVTK,compute_DiceSimilarity,fid_out);
    
    disp('Scaling the transformation field back to its original size...');
    Vx = Vx_new - Vx_ident;
    Vy = Vy_new - Vy_ident;
    Vz = Vz_new - Vz_ident;
    
    % upsample
    Vx = resize3Dmatrix(ori_size_x,ori_size_y,ori_size_z,Vx/scale);
    Vy = resize3Dmatrix(ori_size_x,ori_size_y,ori_size_z,-Vy/scale);
    Vz = resize3Dmatrix(ori_size_x,ori_size_y,ori_size_z,Vz/scale);
    
    I0 = resize3Dmatrix(ori_size_x,ori_size_y,ori_size_z,I0);
   
    disp('Scaling the image back to its original size...');
    count=count+1;
    
end

if compute_DiceSimilarity 
    
    if PlotAll==0
    figure
    nz=size(I1,3);
    position_slice_z=round(nz/2);
    [ I1label ] = ConvertIntensity( I1 );
    subplot(1,2,1)
    imagesc(flipud(I1label(:,:,position_slice_z)'))
    title(['Slice of Fixed Image at z= ', num2str(position_slice_z)])
    axis image
    colormap('gray')
    drawnow
 
    nz=size(I0,3);
    position_slice_z=round(nz/2);
    [ I0label ] = ConvertIntensity( I0 );
    subplot(1,2,2)
    imagesc(flipud(I0label(:,:,position_slice_z)'))
    title(['Slice of registered Moving Image at z= ', num2str(position_slice_z)])
    axis image
    colormap('gray')
    drawnow
    end
    
    fprintf(['Output written to ', output_file])
    fclose(fid_out);
 
end

if SaveVideo
    frameRate = 15;
    vidObj = VideoWriter('RegistrationExample3D.avi');
    vidObj.FrameRate = frameRate;
    open(vidObj)
    writeVideo(vidObj, f)
    close(vidObj)
end
