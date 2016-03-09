% Generates the .mex files using codegen


cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = true;

%% Define argument types for entry-point 'BsplineCompose2D'.
ARGS = cell(1,1);
ARGS{1} = cell(6,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{5} = coder.typeof(0);
ARGS{1}{6} = coder.typeof(0);

%% Invoke MATLAB Coder.
disp('Generating BsplineCompose2D...')
codegen -config cfg BsplineCompose2D -args ARGS{1} -o BsplineCompose2D 



cfg = coder.config('mex');
cfg.GenerateReport = true;

%% Define argument types for entry-point 'BsplineComposeImage2D'.
ARGS = cell(1,1);
ARGS{1} = cell(5,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]);
ARGS{1}{4} = coder.typeof(0);
ARGS{1}{5} = coder.typeof(0);

%% Invoke MATLAB Coder.
disp('Generating BsplineComposeImage2D...')
codegen -config cfg BsplineComposeImage2D -args ARGS{1} -o BsplineComposeImage2D


%% Define argument types for entry-point 'img2coef2D'.
ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]);

%% Invoke MATLAB Coder.
disp('Generating img2coef2D...')
codegen -config cfg img2coef2D -args ARGS{1} -o img2coef2D

