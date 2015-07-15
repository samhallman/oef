function [E,Es] = detect( I, model )
% [E,Es] = detect(I, model)
%
% INPUTS
%  I     - HxWx3 RGB image
%  model - OEF model (output of train.m)
%
% OUTPUTS
%  E     - HxW image of boundary strengths in [0,inf)
%  Es    - HxWxO boundary strengths at different orientations,
%          where O = model.opts.nOrients is set at train time.
%
% Note that you probably want to set model.opts.collapse=1.
% This dramatically speeds up detection at almost no cost in
% accuracy (ODS/AP drop by < .01). You can also get slightly
% better AP by setting model.opts.calibrate = 1, but this
% requires calibration weights model.beta (see calibrate.m).

% check parameters
opts = model.opts; k=length(opts.scales);
assert(length(opts.sharpen) == k);
assert(opts.stride >= opts.shrink);
assert(opts.nTreesEval <= opts.nTrees);
if(opts.calibrate), assert(length(model.beta)==k); end
if(~isfield(opts,'threshW')), model.opts.threshW=0; end

% run ssdetect once for each scale
scales=opts.scales; sharpen=opts.sharpen;
beta=model.beta; [h,w,~]=size(I); Es=0;
for i=1:k, sc=scales(i); sh=sharpen(i);
  if(opts.calibrate), model.beta=beta(i); end
  I1=imResample(I,sc); model.opts.sharpen=sh;
  E1=ssdetect(I1,model); Es=Es+imResample(E1,[h w]);
end; Es=Es/k;

% flatten and optionally perform nms
if ~opts.nms, E=sum(Es,3);
else E=nms(Es,1,80,opts.nThreads); end

end

function E = ssdetect( I, model )

% pad image, making divisible by 4
[h,w,~]=size(I); opts=model.opts;
I = impad4(I,opts.imWidth/2);

% apply forest to image and build weights W
[chnsReg,chnsSim] = edgesChns( I, opts );
if(opts.sharpen), I=convTri(single(I),1); end
ind = imToInd(model,chnsReg,chnsSim);
W = indToWavg(ind,model);
if(opts.threshW>0), W(W<opts.threshW)=0; end
if(opts.calibrate), W=1-exp(-model.beta*W); end
% optionally collapse d~=0 channels onto d=0
if(opts.collapse), V=collapse(W,opts); W(:)=0;
  n=opts.nDists; m=opts.nOrients; cen=(n+1)/2;
  inds=reshape(1:m*n,n,m); W(:,:,inds(cen,:))=V;
end

% build E from distributions W
E = WtoE(model,W,I);
r = opts.gtWidth/2;
E = E(1+r:h+r,1+r:w+r,:);

end
