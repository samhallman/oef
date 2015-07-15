function im = angleimage(bmap,r,show)

imsize = size(bmap);
edgelist = edgelink(bmap);

im = nan(imsize);
for i = 1:length(edgelist)
  x = edgelist{i}(:,2);
  y = edgelist{i}(:,1);
  theta = contour_angles(x,y,r)*180/pi;
  theta(theta<0) = theta(theta<0)+360;
  im( sub2ind(imsize,y,x) ) = theta;
end

if nargin>3 && show
  figure, imagesc(im), colorbar
  colormap([0 0 0;jet])
  title(['r=' num2str(r)])
end
