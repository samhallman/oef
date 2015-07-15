function bmap = cleanbmap( bmap )

list = edgelink(bmap);
list = cleanedgelist(list,7);
bmap = edgelist2image(list,size(bmap));
bmap = bwmorph(bwmorph(bmap,'fill'),'thin',inf);

% [update 17 Sep 2014]
% If the output bmap has a single, isolated 'on' pixel
% then edgelink will not find it, which leads to bugs
bmap = bwmorph(bmap,'clean');
