function [ dxKernel, dyKernel, dzKernel ] = DiffBsplineKernel3D( )
%create 4x4x4 kernels used for evaluate partial derivative in x , y, and Z direction respectively.    

%mid point value at midpoint
b=[1/48, 23/48,  23/48, 1/48];
%diff value
db_x=[-1/8 -5/8, 5/8 1/8];
db_y=[1/8, 5/8, -5/8, -1/8];
db_z=[-1/8,-5/8, 5/8, 1/8];

dxKernel=zeros(4,4,4);
dyKernel=zeros(4,4,4);
dzKernel=zeros(4,4,4);


for i=1:4
    for j=1:4
        for k=1:4
            dxKernel(j,i,k)=db_x(i)*b(j)*b(k);
        end
    end
end
   
for i=1:4
    for j=1:4
        for k=1:4
            dyKernel(j,i,k)=b(i)*db_y(j)*b(k);
        end
    end
end

for i=1:4
    for j=1:4
        for k=1:4
            dzKernel(j,i,k)=b(i)*b(j)*db_z(k);
        end
    end
end

end

