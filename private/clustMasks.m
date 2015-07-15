function [S,E] = clustMasks( opts )

gtWidth = opts.gtWidth;
nDists = opts.nDists;
nOrients = opts.nOrients;
if rem(gtWidth,2) ~= 0,
  error('gtWidth should be even');
end

% define gtWidth x gtWidth grid
[xx,yy] = meshgrid(-gtWidth/2+1:gtWidth/2);

% define the distances and orientations
k = gtWidth/2-1;
dist = linspace(k,-k,nDists);
theta = 0:(180/nOrients):(180-.01);

% render seg masks for each cluster
K = nDists*nOrients;
S = zeros(gtWidth,gtWidth,K);
for i = 1:nOrients
  t = theta(i);
  w = [cosd(t);sind(t)];
  zz = w(1)*xx + w(2)*yy;
  for j = 1:nDists
    k = (i-1)*nDists+j;
    S(:,:,k) = zz > dist(j);
  end
end
% check for bugs
assert(all(sum(sum(S)) < gtWidth^2));

% demonstrate how to convert segs S to edges E
if nargout > 1
  E = zeros(size(S));
  for k = 1:K
    E(:,:,k) = gradientMag(single(S(:,:,k)))>.01;
  end
end
