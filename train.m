function model = train( varargin )
% Train on OEF model
%  opts = train()
%  model = train(opts)
%
% Some important options are
%  nTrees - the number of trees to train
%  nNeg - the number of background patches
%  nPos - the number of non-background patches. Each has a
%    ground truth distance (d) and orientation (theta), and
%    we sample over classes (d,theta) equally
%  modelFnm - a name, e.g. 'model'
%  bsdsDir - the path to the BSDS500 dataset
%  useParfor - if true then train trees in parallel
%
% NOTE that the calibration code contains a number of hard-
% coded constants, and will break if you change some of the
% options not mentioned above (e.g. gtWidth, shrink). 
%
% EXAMPLE
%   opts=train; opts.nPos=1e6; opts.nNeg=1e6; opts.nTrees=12;
%   opts.useParfor=1; opts.bsdsDir='../bsds/'; model=train(opts);

% get default parameters
dfs={'imWidth',32, 'gtWidth',16, 'nOrients',8, 'nDists',15, ...
  'nPos',2e6, 'nNeg',2e6, 'nTrees',24, 'fracFtrs',1/8, ...
  'shrink',2, 'minCount',1, 'minChild',8, 'maxDepth',64, ...
  'split','gini', 'nHistBins',4, 'grdSmooth',0, 'chnSmooth',2, ...
  'simSmooth',8, 'normRad',4, 'nCells',5, 'angleRad',6, ...
  'seed',1, 'calibrate',1, 'useParfor',0, 'cacheDir','cache/', ...
  'modelFnm','model', 'bsdsDir','~/bsds/'};
opts = getPrmDflt(varargin,dfs,1);
if(nargin==0), model=opts; return; end

% if forest exists load it and return
forestDir = [opts.cacheDir '/forest/'];
forestFnm = [forestDir opts.modelFnm];
if(exist([forestFnm '.mat'], 'file'))
  load([forestFnm '.mat']); return; end

% compute constants and store in opts
nTrees=opts.nTrees; nCells=opts.nCells; shrink=opts.shrink;
nOrients=opts.nOrients; nDists=opts.nDists;
angleRad=opts.angleRad; opts.nClusts=nOrients*nDists;
opts.nPos=round(opts.nPos); opts.nNeg=round(opts.nNeg);
imWidth=opts.imWidth; gtWidth=opts.gtWidth;
imWidth=round(max(gtWidth,imWidth)/shrink/2)*shrink*2;
opts.imWidth=imWidth; opts.gtWidth=gtWidth;
nChnsGrad=(opts.nHistBins+1)*2; nChnsColor=3;
nChns = nChnsGrad+nChnsColor; opts.nChns = nChns;
opts.nChnFtrs = imWidth*imWidth*nChns/shrink/shrink;
opts.nSimFtrs = (nCells*nCells)*(nCells*nCells-1)/2*nChns;
opts.nTotFtrs = opts.nChnFtrs + opts.nSimFtrs;
opts.clustFnm = sprintf('clusters_%d_%d_%d_%d_%d_%d', ...
  imWidth, gtWidth, nDists, nOrients, shrink, angleRad); disp(opts);

% compute cluster data if it doesn't exist
clustFnm = [opts.cacheDir '/clust/' opts.clustFnm '.mat'];
if ~exist(clustFnm,'file'), disp('Clust file not found');
  clusters = genClustData(opts); save(clustFnm,'clusters');
end

% train nTrees random trees (can use parfor if enough memory)
stream = RandStream('mrg32k3a','Seed',opts.seed);
if(opts.useParfor), parfor i=1:nTrees, trainTree(opts,stream,i); end
else for i=1:nTrees, trainTree(opts,stream,i); end; end

% build model struct
model = mergeTrees(opts); model.beta = [];
% add test time options (can be changed later)
opts=model.opts; opts.stride=opts.shrink;
opts.nThreads=4; opts.nTreesEval=round(nTrees/2);
opts.scales=[.25 .5 1 2]; opts.sharpen=[1 1 2 2];
opts.collapse=1; opts.nms=0; model.opts=opts;
% learn calibration weights model.beta
if(opts.calibrate), model = calibrate(model); end
% save model
if(~exist(forestDir,'dir')), mkdir(forestDir); end
save([forestFnm '.mat'], 'model', '-v7.3');

end

function trainTree( opts, stream, treeInd )
% Train a single tree in forest model.

% location of ground truth
trnImgDir = [opts.bsdsDir '/images/train/'];
trnGtDir = [opts.bsdsDir '/groundTruth/train/'];
imgIds=dir([trnImgDir '*.jpg']); imgIds={imgIds.name};
nImgs=length(imgIds); for i=1:nImgs, imgIds{i}=imgIds{i}(1:end-4); end

% extract commonly used options
imWidth=opts.imWidth; imRadius=imWidth/2;
gtWidth=opts.gtWidth; gtRadius=gtWidth/2;
nChns=opts.nChns; nTotFtrs=opts.nTotFtrs;
nPos=opts.nPos; nNeg=opts.nNeg; shrink=opts.shrink;

% finalize setup
treeDir = [opts.cacheDir '/tree/'];
treeFn = [treeDir opts.modelFnm '_tree'];
if(exist([treeFn int2str2(treeInd,3) '.mat'],'file')), return; end
fprintf('\n-------------------------------------------\n');
fprintf('Training tree %d of %d\n',treeInd,opts.nTrees); tStart=clock;

% set global stream to stream with given substream (will undo at end)
streamOrig = RandStream.getGlobalStream();
set(stream,'Substream',treeInd);
RandStream.setGlobalStream( stream );

% load the cluster data
clstr = load([opts.cacheDir '/clust/' opts.clustFnm]);
clstr = clstr.clusters; nClusts = size(clstr.clusts,3);

% sample n patch locations per cluster
n = floor(nPos/nClusts); centers = [];
for c = 1:nClusts, ids = find(clstr.clustId==c);
  ids = ids(randperm(length(ids),min(n,length(ids))));
  centers = [centers; [clstr.x(ids) clstr.y(ids) ...
    clstr.imId(ids) clstr.clustId(ids) clstr.gtId(ids)]];
end; clear clstr;

% collect positive and negative patches and compute features
fids=sort(randperm(nTotFtrs,round(nTotFtrs*opts.fracFtrs)));
ftrs = zeros(nPos+nNeg,length(fids),'single');
labels = zeros(nPos+nNeg,1); k = 0;
tid = ticStatus('Collecting data',1,1);
for i = 1:nImgs
  % get image and compute channels
  gt=load([trnGtDir imgIds{i} '.mat']); gt=gt.groundTruth;
  I=imread([trnImgDir imgIds{i} '.jpg']); siz=size(I);
  [chnsReg,chnsSim] = edgesChns(impad4(I),opts);
  % sample positive and negative locations
  nGt=length(gt); xy=[]; k1=0; B=false(siz(1),siz(2));
  B(shrink:shrink:end,shrink:shrink:end)=1;
  B([1:imRadius end-imRadius:end],:)=0;
  B(:,[1:imRadius end-imRadius:end])=0;
  % sample positive locations
  centers1=centers(centers(:,3)==i,:);
  x=centers1(:,1); y=centers1(:,2); k2=length(y);
  xy=[xy; x y centers1(:,4)]; k1=k1+k2;
  % sample negative locations
  for j=1:nGt
    M=gt{j}.Boundaries; M(bwdist(M)<gtRadius)=1;
    [y,x]=find(~M.*B); k2=min(length(y),ceil(nNeg/nImgs/nGt));
    rp=randperm(length(y),k2); y=y(rp); x=x(rp);
    xy=[xy; x y ones(k2,1)*(nClusts+1)]; k1=k1+k2; %#ok<AGROW>
  end
  if(k1>size(ftrs,1)-k), k1=size(ftrs,1)-k; xy=xy(1:k1,:); end
  % crop patches
  psReg=zeros(imWidth/shrink,imWidth/shrink,nChns,k1,'single');
  psSim=psReg; ri=imRadius/shrink;
  for j=1:k1, xy1=xy(j,:); xy2=round(xy1/shrink);
    psReg(:,:,:,j)=chnsReg(xy2(2)-ri+1:xy2(2)+ri,xy2(1)-ri+1:xy2(1)+ri,:);
    psSim(:,:,:,j)=chnsSim(xy2(2)-ri+1:xy2(2)+ri,xy2(1)-ri+1:xy2(1)+ri,:);
  end
  % compute features and store
  ftrs1=[reshape(psReg,[],k1)' stComputeSimFtrs(psSim,opts)];
  ftrs(k+1:k+k1,:)=ftrs1(:,fids); labels(k+1:k+k1)=xy(:,end);
  k=k+k1; if(k==size(ftrs,1)), tocStatus(tid,1); break; end
  tocStatus(tid,i/nImgs);
end
if(k<size(ftrs,1)), ftrs=ftrs(1:k,:); labels=labels(1:k); end

% train structured edge classifier (random decision tree)
pTree=struct('minCount',opts.minCount, 'minChild',opts.minChild, ...
  'maxDepth',opts.maxDepth, 'split',opts.split, 'H',nClusts+1);
tree=forestTrain(ftrs,labels,pTree);
tree.fids(tree.child>0) = fids(tree.fids(tree.child>0)+1)-1;
if(~exist(treeDir,'dir')), mkdir(treeDir); end
save([treeFn int2str2(treeInd,3) '.mat'],'tree'); e=etime(clock,tStart);
fprintf('Training of tree %d complete (time=%.1fs).\n',treeInd,e);
RandStream.setGlobalStream( streamOrig );

end

function ftrs = stComputeSimFtrs( chns, opts )
% Compute self-similarity features (order must be compatible w mex file).
w=opts.imWidth/opts.shrink; n=opts.nCells; if(n==0), ftrs=[]; return; end
nSimFtrs=opts.nSimFtrs; nChns=opts.nChns; m=size(chns,4);
inds=round(w/n/2); inds=round((1:n)*(w+2*inds-1)/(n+1)-inds+1);
chns=reshape(chns(inds,inds,:,:),n*n,nChns,m);
ftrs=zeros(nSimFtrs/nChns,nChns,m,'single');
k=0; for i=1:n*n-1, k1=n*n-i; i1=ones(1,k1)*i;
  ftrs(k+1:k+k1,:,:)=chns(i1,:,:)-chns((1:k1)+i,:,:); k=k+k1; end
ftrs = reshape(ftrs,nSimFtrs,m)';
end
