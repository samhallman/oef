function W = buildW( I, model, method )
% W = buildW( I, model, [method] )
% method is either 'avg' or 'vote'
% This file is used for calibration

if nargin < 3, method = 'avg'; end

% Select method for building W from tree predictions
switch method
  case 'avg', indToW = @indToWavg;
  case 'vote', indToW = @indToWvote;
end
% Run forest over the image and build W
[chns,chnsSs] = edgesChns( impad4(I), model.opts );
W = indToW( imToInd(model,chns,chnsSs), model );
