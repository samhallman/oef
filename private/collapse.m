function V = collapse( W, opts, usemex )
% V = collapse( W, opts, usemex )

if nargin < 3, usemex = 1; end

% TODO get rid of these restrictions
assert(opts.nOrients==8);
assert(opts.stride==opts.shrink);

orthog = [90:-22.5:-67.5]-90;
if usemex, V = collapseMex( W, single(orthog), opts.gtWidth, ...
  opts.nDists, opts.nOrients, opts.shrink, opts.nThreads);
else
  nDists = opts.nDists;
  V = 0*W(:,:,1:opts.nOrients);
  for o = 1:opts.nOrients,
    Wori = W(:,:,(o-1)*nDists+1:o*nDists);
    V(:,:,o) = collapse_orient( Wori, orthog(o), opts );
  end
end

function W = collapse_orient( Ws, theta, opts )
nDists=opts.nDists; shrink=opts.shrink;
% pad chns and define central pixel grid 
k=(nDists-1)/2; cen=k+1; r=opts.gtWidth/2-1;
p=ceil(r/shrink); Ws=imPad(Ws,p,'replicate');
[h,w,~]=size(Ws); xx=1+p:w-p; yy=1+p:h-p;
[X,Y]=meshgrid(xx,yy); W=Ws(yy,xx,cen);
% make step vector in orthog direction
dx=cosd(theta); dy=sind(theta);
m=max(abs([dx dy])); dx=dx/m; dy=dy/m;
% translate d~=0 channels to line up with d=0
dists = linspace(-r,r,nDists);
for i=1:nDists, d=dists(i); if(d==0), continue; end
  Xq = X+dx*d/shrink; Yq = Y-dy*d/shrink;
  W = W + interpW_faster( Ws(:,:,i), Xq, Yq );
end

function W = interpW( W, x, y )
x0=floor(x); x1=x0+1; dx0=x-x0; dx1=1-dx0;
y0=floor(y); y1=y0+1; dy0=y-y0; dy1=1-dy0; h=size(W,1);
W = W((x0-1)*h+y0).*dx1.*dy1 + W((x1-1)*h+y0).*dx0.*dy1 + ...
    W((x0-1)*h+y1).*dx1.*dy0 + W((x1-1)*h+y1).*dx0.*dy0;

function W = interpW_faster( W, x, y )
% avoid expensive interp math for integer input
x0 = floor(x); intx = all(all(x==x0));
y0 = floor(y); inty = all(all(y==y0)); h = size(W,1);
if( intx && inty ), W = W((x-1)*h+y);
elseif( intx ), y1=y0+1; dy0=y-y0; dy1=1-dy0;
  W = W((x-1)*h+y0).*dy1 + W((x-1)*h+y1).*dy0;
elseif( inty ), x1=x0+1; dx0=x-x0; dx1=1-dx0;
  W = W((x0-1)*h+y).*dx1 + W((x1-1)*h+y).*dx0;
else x1=x0+1; dx0=x-x0; dx1=1-dx0; y1=y0+1; dy0=y-y0; dy1=1-dy0;
  W = W((x0-1)*h+y0).*dx1.*dy1 + W((x1-1)*h+y0).*dx0.*dy1 + ...
      W((x0-1)*h+y1).*dx1.*dy0 + W((x1-1)*h+y1).*dx0.*dy0;
end

