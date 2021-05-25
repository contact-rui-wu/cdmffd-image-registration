function I0 = applyTransform2Image2D(I0,Vx,Vy)

% I00 = I0;
% scale = 2^-(level-1); % level=1 if we're just applying the transform
% I0 = imresize(I0,scale);
% I0_ori=imresize(I00,scale);
% (in original code, images are scaled depending on level; here level=1 so
% no need scaling)
I0_ori = I0;

%-----begin: regImage2D-----

ny = size(I0,1);
nx = size(I0,2);

II0_ori = I0_ori(:)';

coef_I0 = img2coef2D(nx,ny,II0_ori);
CI0 = coef_I0(:)';

[Vx_ident, Vy_ident] = meshgrid((0.5:nx-0.5), (0.5:ny-0.5));

% Vx = Vx_new - Vx_ident;
% Vy = Vy_new - Vy_ident;
% Vy = -Vy;
% (in original code, Vx_new and Vy_new are computed by BsplineCompose2D)
Vx_new = Vx + Vx_ident;
Vy_new = -Vy + Vy_ident;

Vx_new = reshape(Vx_new, 1, nx*ny);
Vy_new = reshape(Vy_new, 1, nx*ny);

I0 = BsplineComposeImage2D(Vx_new, Vy_new, CI0, nx, ny);

%-----end: regImage2D-----

% % scaling the transformation field back to its original size
% % upsample
% Vx = imresize(Vx/scale,size(I00));
% Vy = imresize(Vy/scale,size(I00));
% % scaling the image back to its original size
% I0 = imresize(I0/scale,size(I00));
% (in original code, images are scaled depending on level; here level=1 so
% no need scaling)

end