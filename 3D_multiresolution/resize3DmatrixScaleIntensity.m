function [I ] = resize3DmatrixScaleIntensity( nx,ny,nz,I )
% resize3Dmatrix  resize 3 dimensional image (I) to a given size
% (nx,ny,nz), and scale the intensity to 0-255 
% Input: I = 3 dimensional matrix. 
%        nx, ny,nz = desired size in x, y ,and z direction.
% Output: I = rezised 3 dimensional matrix. 

[y, x, z]=...
   ndgrid(linspace(1,size(I,1),ny),...
          linspace(1,size(I,2),nx),...
          linspace(1,size(I,3),nz));

I=interp3(I,x,y,z);
I=(I-min(I(:)))/(max(I(:))-min(I(:)))*255;
I=round(I);

end

