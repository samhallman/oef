function orient = binangles360(angles,nBins)
% Assume angles lie in [0,360]

if any(angles<0 | angles>360)
  error('angles should be in the range [0,360]');
end

%
% change from [0,360] to [-270,90]
%
angles(angles>90) = angles(angles>90)-360;

delta = 360/nBins;
a = -270 + delta/2;
b =  +90 - delta/2;
edges = a:delta:b;
[~,orient] = histc(angles,edges);
orient = nBins-orient+1;
orient(orient == nBins+1) = 1;
