function [ Vx_new, Vy_new, Vz_new ] = BsplineCompose3D( nx,ny,nz,VX, VY, VZ, CX, CY, CZ)
%compose the transformation field Vx, Vy, Vz with the update field given by in terms of coefficient (CX, CY, and CZ).
%Note: Vx, Vy Vy should be between 0...nx, 0..ny and 0..nz respectively

%input: nx,ny,nz = size of image in x, y, and z direction respectively
%       VX,VY,VZ = transformation field in x, y, and z direction respectively
%       CX, CY, CZ = coefficient of transformation field in  x, y, and z direction respectively

%output: Vx_new, Vy_new, Vz_new = composed transformation field in x, y, and z direction respectively

%degree of B-spline level set function
p=3;
q=3;
r=3;

Vx=reshape(VX,[ny nx nz]);
Vy=reshape(VY,[ny nx nz]);
Vz=reshape(VZ,[ny nx nz]);

coef_x=reshape(CX,[ny+q nx+p nz+r]);
coef_y=reshape(CY,[ny+q nx+p nz+r]);
coef_z=reshape(CZ,[ny+q nx+p nz+r]);

Vx_new = zeros(size(Vx));
Vy_new = zeros(size(Vy));
Vz_new = zeros(size(Vz));
[Vx_ident, Vy_ident, Vz_ident] = meshgrid((0.5:nx-0.5),(0.5:ny-0.5),(0.5:nz-0.5));
parfor k=1:nz
    for j=1:ny
        for i=1:nx
            
         
            x_disp = Vx(j,i,k);
            y_disp = Vy(j,i,k);
            z_disp = Vz(j,i,k);
            %keep the displacement within the bounds of the image
            x_disp = max(0, x_disp);
            y_disp = max(0, y_disp);
            z_disp = max(0, z_disp);
            
            x_disp = min(x_disp, nx-eps);
            y_disp = min(y_disp, ny-eps);
            z_disp = min(z_disp, nz-eps);
            
            k1 = ceil(x_disp-2+eps);
            l1 = ceil(y_disp-2+eps);
            m1 = ceil(z_disp-2+eps);
            
            x = x_disp-k1;
            y = y_disp-l1;
            z = z_disp-m1;
            
            basis_x = [(2-x)^3/6, 2/3-(x-1)^2+(x-1)^3/2, 2/3-(x-2)^2+(2-x)^3/2, (x-1)^3/6];
            basis_y = [(2-y)^3/6, 2/3-(y-1)^2+(y-1)^3/2, 2/3-(y-2)^2+(2-y)^3/2, (y-1)^3/6];
            basis_z = [(2-z)^3/6, 2/3-(z-1)^2+(z-1)^3/2, 2/3-(z-2)^2+(2-z)^3/2, (z-1)^3/6];
            
            %basis_y = basis_y(end:-1:1);
            coef_loc_x = coef_x(l1+2:l1+5,k1+2:k1+5,m1+2:m1+5);
            coef_loc_y = coef_y(l1+2:l1+5,k1+2:k1+5,m1+2:m1+5);
            coef_loc_z = coef_z(l1+2:l1+5,k1+2:k1+5,m1+2:m1+5);

            Vx_temp = 0;
            for u=1:4
                for v=1:4
                    for w=1:4
                        Vx_temp=Vx_temp+coef_loc_x(v,u,w)*basis_x(u)*basis_y(v)*basis_z(w);                       
                    end
                end
            end
            Vx_new(j,i,k)=Vx_temp;
            
            Vy_temp = 0;
            for u=1:4
                for v=1:4
                    for w=1:4
                        Vy_temp=Vy_temp + coef_loc_y(v,u,w)*basis_x(u)*basis_y(v)*basis_z(w);
                        
                    end
                end
            end
            Vy_new(j,i,k)=Vy_temp;
            
            
            Vz_temp = 0;
            for u=1:4
                for v=1:4
                    for w=1:4
                        Vz_temp=Vz_temp+coef_loc_z(v,u,w)*basis_x(u)*basis_y(v)*basis_z(w);
                        
                    end
                end
            end
            Vz_new(j,i,k)=Vz_temp;
            
        end
    end
end

Vx_new(:,1:3,:)=Vx_ident(:,1:3,:);
Vx_new(:,nx-2:nx,:)=Vx_ident(:,nx-2:nx,:);
Vy_new(1:3,:,:)=Vy_ident(1:3,:,:);
Vy_new(ny-2:ny,:,:)=Vy_ident(ny-2:ny,:,:);
Vz_new(:,:,1:3)=Vz_ident(:,:,1:3);
Vz_new(:,:,nz-2:nz)=Vz_ident(:,:,nz-2:nz);