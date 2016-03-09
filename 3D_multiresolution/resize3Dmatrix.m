function [I ] = resize3Dmatrix( nx,ny,nz,I )
% resize3Dmatrix resize 3 dimensional matrix (I) to a given size (nx,ny,nz)
% Input: I = 3 dimensional matrix. 
%        nx, ny,nz = desired size in x, y ,and z direction.
% Output: I = rezised 3 dimensional matrix. 

[y, x, z]=...
   ndgrid(linspace(1,size(I,1),ny),...
          linspace(1,size(I,2),nx),...
          linspace(1,size(I,3),nz));
I=interp3(I,x,y,z);
end

