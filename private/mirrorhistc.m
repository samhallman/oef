function [n,b] = mirrorhistc( x, e )
% [N,B] = mirrorhistc( X, E )
%
% Suppose E = [a b]. Then the binning
% used is (-b,-a] (-a,a) [a,b).
% 
% Suppose E = [a b c]. Then the binning
% used is (-c,-b] (-b,-a] (-a,a) [a,b) [b,c).
%
% ...and so on.

if any(e<=0) || ~all(diff(e)>0)
  error('invalid bin edges')
end

if abs(max(x)) >= e(end)
  error('data is out of bounds');
end

if ndims(x) > 2 || min(size(x)) > 1
  error('data should be a vector');
end

x = x(:)';
e = e(:)';

e = [0 e];
lastBin = length(e)-1;

[nn,bn] = histc(-x(x<0),e);
bn = lastBin-bn+1;
nn = nn(end-1:-1:1);

[np,bp] = histc(x(x>=0),e);
bp = bp+lastBin-1;
np = np(1:end-1);

b = 0*x; b(x<0) = bn; b(x>=0) = bp;
n = [nn(1:end-1) nn(end)+np(1) np(2:end)];
