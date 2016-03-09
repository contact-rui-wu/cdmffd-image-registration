% Generates the .mex files using codegen


cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = true;

%% Define argument types for entry-point 'BsplineCompose3D'.
ARGS = cell(1,1);
ARGS{1} = cell(9,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{6} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{7} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{8} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{9} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
disp('Generating BsplineCompose3D...')
codegen -config cfg BsplineCompose3D -args ARGS{1} -o BsplineCompose3D 



cfg = coder.config('mex');
cfg.GenerateReport = true;

%% Define argument types for entry-point 'BsplineComposeImage3D'.
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{6} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{7} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
disp('Generating BsplineComposeImage3D...')
codegen -config cfg BsplineComposeImage3D -args ARGS{1} -o BsplineComposeImage3D


%% Define argument types for entry-point 'img2coef3D'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0);

%% Invoke MATLAB Coder.
disp('Generating img2coef3D...')
codegen -config cfg img2coef3D -args ARGS{1} -o img2coef3D

