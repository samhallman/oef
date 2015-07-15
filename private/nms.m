function E = nms( Es, alpha, tol, nThreads )
% E = nms( Es, alpha, tol, nThreads )

if(nargin<2), alpha=1; end
if(nargin<3), tol=80; end
if(nargin<4), nThreads=4; end
nOrients=size(Es,3); tol=tol*pi/180;

% Flatten and smooth slightly (it helps)
E = sum(Es,3);
E = convTri(E,1);

% Compute orientation map
[~,Ot] = max(Es,[],3);
Ot = single(Ot-1)/nOrients*pi;
Oe = edgeOrient(E,4);
if(tol>0), outlier=abs(Oe-Ot)>tol;
  Oe(outlier)=Ot(outlier); end
Xt = cos(Ot); Yt = sin(Ot);
Xe = cos(Oe); Ye = sin(Oe);
X = Xt + alpha*(Xe-Xt);
Y = Yt + alpha*(Ye-Yt);
O = atan( Y ./ X );

% Use those orientations to do nms
E = edgesNmsMex(E,O,1,5,1.01,nThreads);

function O = edgeOrient( W, r )
% compute approximate orientation map
[Ox,Oy] = gradient2(convTri(W,r));
[Oxx,~] = gradient2(Ox); [Oxy,Oyy] = gradient2(Oy);
O = mod( atan( Oyy.*sign(-Oxy)./(Oxx+1e-5) ), pi );
