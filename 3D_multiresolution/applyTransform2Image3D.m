function I0 = applyTransform2Image3D(I0,Vx,Vy,Vz)

ny = size(I0,1);
nx = size(I0,2);
nz = size(I0,3);

[Vx_ident,Vy_ident,Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));

coef_I0 = img2coef3D(I0(:)',nx,ny,nz);

Vx_new = Vx + Vx_ident;
Vy_new = -Vy + Vy_ident;
Vz_new = Vz + Vz_ident;

VX_new = Vx_new(:)';
VY_new = Vy_new(:)';
VZ_new = Vz_new(:)';

CI0 = coef_I0(:)';

I0 = BsplineComposeImage3D(nx,ny,nz,VX_new,VY_new,VZ_new,CI0);

end