function angles = contour_angles(x,y,r)

% edgelink can return one-pixel "edges"
if length(x) == 1, angles = 0;
  warning('edge is length 1'); return;
end

if length(x) < 3
  dy = y(end)-y(1); dx = x(end)-x(1);
  assert( dy~=0 || dx~=0 );
  angles = atan2(-dy,dx)*ones(size(x));
elseif length(x) < 2*r+1
  angles = contour_angles(x,y,r-1);
else
  if x(1)==x(end) && y(1)==y(end)
    angles = loop_angles(x,y,r);
  else
    angles = nonloop_angles_fast(x,y,r);
    angles(1:r) = angles(r+1);
    angles(end-r+1:end) = angles(end-r);
  end
end


function theta = loop_angles(x,y,r)

l = 2*r+1;
t = -floor(l/2):floor(l/2);
T = [ones(size(t)); t; t.^2]';
T = pinv(T);

n = length(x);
theta = zeros(n,1);
for c = 1:n
  wind = c-r:c+r;
  wind = 1+mod(wind-1,n);
  % fit polynomials around x(c),y(c)
  a = T*x(wind);
  b = T*y(wind);
  theta(c) = atan2(a(2),b(2))+(pi/2);
end


function theta = nonloop_angles_fast(x,y,r)

% Handle the non-border cases quickly with convolution
t = [-r:r]';
T = [ones(size(t)) t t.^2];
T = pinv(T); kern = T(2,:);
a = conv(x,kern,'same');
b = conv(y,kern,'same');
theta = atan2(a,b)+(pi/2);
