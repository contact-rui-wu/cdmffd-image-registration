function [I] = readrawb(filename)
%reads the binaries from BrainWeb (raw byte format, not raw short!)

x_size = 362;
y_size = 434;
z_size = 362;

%read the file as 8-bit unsigned integers (range 0...255)
fid = fopen(filename, 'r'); 
I = fread(fid, inf, 'uint8');
I = reshape(I,x_size,y_size,z_size);

