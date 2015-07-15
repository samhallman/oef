% Quick OEF demo script
% =====================

% Load model
file = 'cache/forest/modelCvpr.mat';
try
  load(file,'model');
catch
  % download the trained model (98 MB)
  unix(['wget --directory-prefix=cache/forest/ ' ...
    'http://www.ics.uci.edu/~shallman/oef/modelCvpr.mat']);
  load(file,'model');
end

% Set some options before calling detect.m
model.opts.nms       = 0;
model.opts.nThreads  = 8;      % set to num logical cpu cores
model.opts.calibrate = true;   % see section 3 of cvpr paper (p.4)
model.opts.collapse  = true;   % see section 5.3 of cvpr paper (pp.6-7)

% Detect boundaries
I = imread('peppers.png');
[E,Es] = detect(I,model);

% View results
figure, subplot(131), imshow(I),
subplot(132), imagesc(E),
colormap gray, axis off image,
subplot(133), nOrients = size(Es,3);
labels = strsplit(num2str(1:nOrients));
montage2( Es, struct('labels',{labels}) );
