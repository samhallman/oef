function [chnsReg,chnsSim] = edgesChns( I, opts )
% Compute feature channels
% [chnsReg,chnsSim] = edgesChns( I, opts )
%
% This is a slightly modified version of the edgesChns
% file from Piotr Dollar's Structured Edge Detection
% Toolbox, https://github.com/pdollar/edges.

shrink=opts.shrink; chns=cell(1,opts.nChns); k=0;
if(size(I,3)==1), cs='gray'; else cs='luv'; end; I=rgbConvert(I,cs);
Ishrink=imResample(I,1/shrink); k=k+1; chns{k}=Ishrink;
for i = 1:2, s=2^(i-1);
  if(s==shrink), I1=Ishrink; else I1=imResample(I,1/s); end
  I1 = convTri( I1, opts.grdSmooth );
  [M,O] = gradientMag( I1, 0, opts.normRad, .01 );
  H = gradientHist( M, O, max(1,shrink/s), opts.nHistBins, 0 );
  k=k+1; chns{k}=imResample(M,s/shrink);
  k=k+1; chns{k}=imResample(H,max(1,s/shrink));
end
chns=cat(3,chns{1:k}); assert(size(chns,3)==opts.nChns);
chnSm=opts.chnSmooth/shrink; if(chnSm>1), chnSm=round(chnSm); end
simSm=opts.simSmooth/shrink; if(simSm>1), simSm=round(simSm); end
chnsReg=convTri(chns,chnSm); chnsSim=convTri(chns,simSm);

end
