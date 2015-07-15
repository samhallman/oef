function N = ncc(seg,width)
% N = ncc( seg, w )
%
% N(i,j) = number of connected components in the
%   w-by-w window centered on location (i,j) of seg
%
% What do I mean by "number of connected components"?  
% Here are some example segmentation patches:
% 
%   1 1 2 2 1 1
%   1 1 2 2 1 1      this patch has THREE segments, even
%   1 1 2 2 1 1      though unique() only returns [1;2].
%   1 1 2 2 1 1
% 
%   3 3 7 7 7 7
%   3 3 7 7 7 7      this patch has FOUR segments, even
%   7 7 3 3 3 3      though unique() is simply [3;7].  
%   7 7 3 3 3 3
%
% THIS CODE IS VERY SLOW, hence all results are cached.

cacheDir = 'cache/clust/nccCache/';
hash = hashMat( [seg(:);width] );
try
  load([cacheDir hash '.mat'],'N');
catch
  siz = size(seg);
  o1 = floor((width-1)/2);
  o2 = ceil((width-1)/2);

  % This takes ~80 seconds per image
  N = nan(siz); B = false(siz);
  B(o1+1:end-o2,o1+1:end-o2) = 1;
  Nu = uniqfilt(seg,o2);
  N(B & Nu==1) = 1;
  [ys,xs] = find(B & Nu>1);
  for i = 1:length(ys)
    y = ys(i); x = xs(i);
    S = seg(y-o1:y+o2,x-o1:x+o2);
    uniq = unique(S); N(y,x) = 0;
    for s = uniq(:)'
      u = unique(bwlabel(S==s,4));
      N(y,x) = N(y,x) + sum(u~=0);
    end
  end

  % convert to uint8 to save disk space
  % (but this converts NaNs to 0s..)
  N = uint8(N);

  save([cacheDir hash '.mat'],'N');
end


function nuniq = uniqfilt(seg,k)
% nuniq = uniqfilt(seg,k)
%
% INPUTS
%   seg     groundTruth{j}.Segmentation
%   k       filter half width
%
% ASSUME THAT SEG LABELS ARE {1,2,...,NSEGS}
%
% This function is a fast version of the following:
% [h w] = size(seg);
% nuniq = zeros(h,w);
% for y = 1:h
%   y1 = max(1,y-k);
%   y2 = min(h,y+k);
%   for x = 1:w
%     x1 = max(1,x-k);
%     x2 = min(w,x+k);
%     patch = seg(y1:y2,x1:x2);
%     nuniq(y,x) = numel(unique(patch));
%   end
% end

% How long does it take to run unique(seg)? I tried this for all
% groundTruth{j}.Segmentation for all images. It took 0.0123 sec
% at most, 0.0027 sec on average, and 0.0022 sec at best.

se = ones(2*k+1);
nuniq = zeros(size(seg));
for i = 1:max(seg(:))
  nuniq = nuniq + imdilate(seg==i,se);
end
