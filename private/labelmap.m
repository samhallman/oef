function [Labels,Dist,Theta] = labelmap( bmap, mask, opts )
% [Labels,Dist,Theta] = labelmap( bmap, mask, opts )
% It is assumed that bmap was already cleaned (cleanbmap.m)

gtWidth = opts.gtWidth;
nOrients = opts.nOrients;
nDists = opts.nDists;
nLabels = nOrients*nDists;
assert(nDists <= gtWidth-1);

% calculate distances and angles
[h,w] = size(bmap);
[Dist,Nearest] = bwdist(bmap);
Nearest = double(Nearest);
Theta = angleimage(bmap,opts.angleRad);
[R,C] = ind2sub([h w],Nearest);
[J,I] = meshgrid(1:w,1:h);
Dx = C-J; Dy = -(R-I);

% focus on locations indicated by mask
locs = find(mask);
dx = Dx(locs); dy = Dy(locs);
nearest = Nearest(locs);
theta = Theta(nearest);
thetaBin = binangles360(theta,nOrients*2);

T = [sind(theta) cosd(theta) 0*theta];
V = [dy dx 0*dx]; U = cross(T,V);
dist = Dist(locs) .* sign(-U(:,3));

% standardize
ind = thetaBin > nOrients;
dist(ind) = -dist(ind);
thetaBin(ind) = thetaBin(ind)-nOrients;

% bin distance
k = (nDists-1)/2;
step = (gtWidth-2)/(nDists-1);
edges = [step*[1:k]-step/2, gtWidth/2];
[~,distBin] = mirrorhistc(dist,edges);

labels = (thetaBin-1)*nDists + distBin';
% IMPORTANT assertion that has caught MANY bugs
assert( all(1<=labels & labels<=nLabels) );

% write result images
Labels=repmat(nLabels+1,[h w]); Dist=Labels; Theta=Labels;
Labels(locs)=labels; Dist(locs)=distBin; Theta(locs)=thetaBin;
