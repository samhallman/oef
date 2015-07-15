function clusters = genClustData( varargin )
% clusters = genClustData( opts )
% Label all 2-segment patches in the train set.
%
% INPUTS
%   imWidth [32], gtWidth [16], nOrients [8],
%   nDists [15], shrink [2], angleRad [6], bsdsDir
%
% OUTPUTS
%   clusters struct w/ fields:
%     x, y, gtId, imId, clustId, clusts, opts

dfs={'imWidth',32, 'gtWidth',16, 'nOrients',8, 'nDists',15, ...
  'shrink',2, 'angleRad',6, 'bsdsDir','~/bsds/'};
opts = getPrmDflt(varargin,dfs,0);

% set up output struct
clusters = struct('x',[], 'y',[], 'gtId',[], ...
  'imId',[], 'clustId',[], 'clusts',[], 'opts',opts);

% set up BSDS train set paths and image ids
imDir = [opts.bsdsDir '/images/train/'];
gtDir = [opts.bsdsDir '/groundTruth/train/'];
ids=dir([imDir '*.jpg']); nImgs=length(ids);
ids={ids.name}; for i=1:nImgs, ids{i}=ids{i}(1:end-4); end

% add path to Peter Kovesi's edge linking code
addpath private/peterkovesi/ % used by labelmap and cleanbmap

% initialize arrays used to compute clust means
gtWidth = opts.gtWidth; gtRadius = gtWidth/2;
nLabels = opts.nOrients * opts.nDists;
counts = zeros(1,nLabels,'uint32');
clusts = zeros(gtWidth,gtWidth,nLabels,'uint32');

tid = ticStatus('Labeling training patches',1,1);

for imgNum = 1:nImgs

  gt = load([gtDir ids{imgNum} '.mat']);
  gt = gt.groundTruth;

  % ignore borders since we can't get patches there,
  % and do not sample points which map to the same
  % locations when images are downsampled by shrink
  locs = false(size(gt{1}.Boundaries));
  s=opts.shrink; r=opts.imWidth/2;
  locs(s:s:end,s:s:end) = 1;
  locs([1:r end-r:end],:) = 0;
  locs(:,[1:r end-r:end]) = 0;
  
  for gtIndex = 1:length(gt)

    seg = gt{gtIndex}.Segmentation;
    bmap = gt{gtIndex}.Boundaries;

    % only label patches with two segments
    bmap = cleanbmap(bmap);
    nsegs = ncc(seg,gtWidth);
    posmask = bwdist(bmap) < gtRadius;
    mask = locs .* posmask .* (nsegs==2);

    % generate output labels
    labels = labelmap(bmap,mask,opts);
    labels = labels(logical(mask));

    [i,j] = find( mask ); nPos = length(i);
    patches = zeros(gtWidth,gtWidth,nPos,'uint8');
    for p = 1:nPos
      is = i(p)-gtRadius+1:i(p)+gtRadius;
      js = j(p)-gtRadius+1:j(p)+gtRadius;
      patches(:,:,p) = bmap(is,js);
    end

    % add these results to the clusters struct
    clusters.x = [clusters.x; j];
    clusters.y = [clusters.y; i];
    clusters.gtId = [clusters.gtId; gtIndex*ones(size(i))];
    clusters.imId = [clusters.imId; imgNum*ones(size(i))];
    clusters.clustId = [clusters.clustId; labels];

    for c = 1:nLabels
      counts(c) = counts(c) + sum(labels==c);
      clusts(:,:,c) = clusts(:,:,c) + ...
        uint32( sum(patches(:,:,labels==c),3) );
    end

  end
  tocStatus(tid,imgNum/nImgs);
end

% divide sums by counts to get means
clusts = double(clusts);
for c = 1:nLabels
  clusts(:,:,c) = clusts(:,:,c) / double(counts(c));
end

clusters.clusts = clusts;
if(0), save('clusters.mat','clusters'); end
