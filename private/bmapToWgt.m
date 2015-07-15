function Wgt = bmapToWgt( bmap, type )
% Wgt = bmapToWgt( bmap, type )
%
% Input bmap should not be padded, and type
% is either 'exact' or 'max' or 'interp'. The
% Wgt returned is binary for 'exact' and 'max'.
%
% The 'max' type seems to work best but has the
% disadvantage that the weights don't sum to 1.
%
% The size of Wgt is the same as the size of
% the W returned by indToW when "ind" was built
% from the original image padded only to have
% height/width divisible by 4 (see imgDemo.m).

bmap = cleanbmap(bmap);
posmask = bwdist(bmap)<(16/2);
L = labelmap( bmap, posmask, struct(...
  'gtWidth',16, 'nDists',15, 'nOrients',8, 'angleRad',6) );

% All of the index math assumes the original image
% was padded to have height/width divisible by 4,
% since that is required by edgesChns
[h,w] = size(bmap);
pad = mod(4-mod([h w],4),4);
h = h+pad(1); w = w+pad(2);

switch type
  % Nearest nbor downsampling
  case 'exact'
    L = L(16:2:h-18,16:2:w-18);
    Wgt = LtoW(L);
  % Charless's idea. Works really well
  case 'max'
    W = LtoW(L);
    is = 16:2:h-18; js = 16:2:w-18;
    Wa = W(is,js,:); Wb = W(is,js+1,:);
    Wc = W(is+1,js,:); Wd = W(is+1,js+1,:);
    Wgt = max(Wa,max(Wb,max(Wc,Wd)));
  % Probably shouldn't use this
  case 'interp'
    Wgt = LtoW(L);
    Wgt = Wgt(16:h-18,16:w-18,:);
    Wgt = imResample(Wgt,0.5,'bilinear');
end

h1 = h/2-16; w1 = w/2-16;
asserteq( size(Wgt), [h1 w1 121] );
