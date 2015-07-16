function [ftrs,labels] = genCalibData( model, varargin )
% [ftrs,labels] = genCalibData( model, opts )
% Train data for model k is X=ftrs(:,:,k), y=labels(:,k)

dfs={'set','train', 'method','avg', 'type','max', ...
  'frac',0.1, 'scale',1, 'radius',0, 'seed',0, 'bsdsDir','~/bsds/'};
opts = getPrmDflt(varargin,dfs,1);
if(nargin==0), ftrs=opts; return; end; disp(opts);

% the Wgt code assumes shrink is 2 (will fix someday)
model.opts.stride=2; assert(model.opts.shrink==2);

% load image ids
imDir = [opts.bsdsDir '/images/' opts.set '/'];
ids=dir([imDir '*.jpg']); ids={ids.name};
nImgs=length(ids); for i=1:nImgs, ids{i}=ids{i}(1:end-4); end

% add path to Peter Kovesi's edge linking code
addpath private/peterkovesi/

% seed the RNG for reproducibility
rng(opts.seed);

parfor i = 1:nImgs
  % load image, resize, and build W
  I = imread([imDir ids{i} '.jpg']);
  I = imResample(I,opts.scale);
  W = buildW( I, model, opts.method );
  
  % generate the ground truth W
  Wgt=idToWgt(ids{i},opts);
  Wgt(:,:,121)=[]; Wgt=single(Wgt);
  [h,w,~]=size(W); Wgt=Wgt(1:h,1:w,:);
  asserteq(size(W),size(Wgt));
  asserteq(size(W,3),size(Wgt,3),120);

  % sample positive and negative locations
  K = 120; P = sum( Wgt(:,:,1:K), 3 );
  [yn,xn] = find(P==0); nNeg = length(yn);
  nPos = numel(P)-nNeg; [h,w,~] = size(W);
  xp = []; yp = []; N = min(nPos,nNeg);
  N = round(opts.frac*N);
  for k = 1:K
    [y,x] = find( Wgt(:,:,k) > 0 );
    n = min(length(y),ceil(N/K));
    rp = randperm(length(y),n);
    yp = [yp; y(rp)]; xp = [xp; x(rp)];
  end
  rp=randperm(nNeg,N); yn=yn(rp); xn=xn(rp);
  ipos = yp+(xp-1)*h; ineg = yn+(xn-1)*h;

  % binarize the problem, and optionally "soften"
  Wgt = single( Wgt > 0 );
  Wgt = convTri( Wgt, opts.radius );
  Wgt = Wgt / max(Wgt(:));

  % collect features and their probabilities
  Wgt = reshape(Wgt, h*w, K);
  W = reshape(W, h*w, K);
  ftrs{i} = W([ipos;ineg],:,:);
  labels{i} = Wgt([ipos;ineg],:);
end

ftrs = cat(1,ftrs{:});
labels = cat(1,labels{:});
