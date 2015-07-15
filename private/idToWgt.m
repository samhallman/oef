function Wgt = idToWgt( id, opts )
% Private function called by genCalibData.
% opts should include: bsdsDir, scale, type

gtDir = [opts.bsdsDir '/groundTruth/all/'];
G = load([gtDir id '.mat']); G = G.groundTruth;
if(opts.scale ~= 1), G = scaleGT( G, opts.scale ); end

% Sum over annotations
Wgt = 0;
for j = 1:length(G)
  Wgt = Wgt + bmapToWgt( G{j}.Boundaries, opts.type );
end

% Normalize!
Wgt = Wgt / length(G);
