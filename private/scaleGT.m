function gt = scaleGT( gt, scale )
% gt = scaleGT( gt, scale )
% Scale the ground truth

for i = 1:length(gt)
  S = gt{i}.Segmentation;
  B = gt{i}.Boundaries;
  % Derive the scaled segmentation
  s = imresize(S,scale,'nearest');
  [~,~,t] = unique(s);
  s = reshape(t,size(s));
  % Derive the corresponding bmap
  b = seg2bdry(s,'imageSize');
  b = bwmorph(b,'thin',inf);
  % Write to struct
  gt{i}.Segmentation = s;
  gt{i}.Boundaries = b;
end
