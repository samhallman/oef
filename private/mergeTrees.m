function model = mergeTrees( opts )
% Merge trees into overall forest model

% load trees from disk
nTrees=opts.nTrees; nNodes=0; nLeaves=0;
treeFn = [opts.cacheDir '/tree/' opts.modelFnm '_tree'];
for i = 1:nTrees
  t=load([treeFn int2str2(i,3) '.mat'],'tree'); t=t.tree;
  if(i==1), trees=t(ones(1,nTrees)); else trees(i)=t; end
  nNodes = max( nNodes, size(trees(i).fids,1) );
  nLeaves = max( nLeaves, sum(trees(i).child==0) );
end

% derive masks for each cluster
[clustSegs,clustEdges] = clustMasks(opts);
nClusts = size(clustEdges,3);

% merge all fields of all trees
model.opts=opts; Z=zeros(nNodes,nTrees,'uint32');
model.thrs=single(Z); model.fids=Z; model.child=Z;
model.hs=Z; model.count=Z; model.depth=Z; model.leafid=Z;
model.distr=zeros(nClusts,nLeaves,nTrees,'single');
for t=1:nTrees, tree=trees(t); nNodes1=size(tree.fids,1);
  model.fids(1:nNodes1,t) = tree.fids;
  model.thrs(1:nNodes1,t) = tree.thrs;
  model.child(1:nNodes1,t) = tree.child;
  model.hs(1:nNodes1,t) = tree.hs;
  model.count(1:nNodes1,t) = tree.count;
  model.depth(1:nNodes1,t) = tree.depth;
  % construct distr without background probabilities
  leaves=tree.child==0; nLeaves1=sum(leaves);
  model.leafid(leaves,t) = 1:nLeaves1;
  model.distr(:,1:nLeaves1,t) = tree.distr(leaves,1:nClusts)';
end
model.hs = uint16(model.hs);
model.segs = uint8(clustSegs);

% store compact representations of sparse binary edge patches
sharpen = 2;
nBnds = sharpen+1;
eBins = cell(nClusts,nBnds);
eBnds = zeros(nClusts,nBnds);
for i = 1:nClusts
  E = clustEdges(:,:,i);
  E0 = 0;
  for j = 1:nBnds
    eBins{i,j} = uint16(find(E & ~E0)'-1);
    E0 = E;
    eBnds(i,j) = length(eBins{i,j});
    E = convTri(single(E),1)>.01;
  end
end
eBins = eBins'; model.eBins = [eBins{:}]';
eBnds = eBnds'; model.eBnds = uint32([0; cumsum(eBnds(:))]);
