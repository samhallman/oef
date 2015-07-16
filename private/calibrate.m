function model = calibrate( model )
% model = calibrate( model )
% Set calibration weights model.beta
%
% This file calls a function that uses a parfor. With
% 12 matlab workers, calibration should take ~20 min.

opts = genCalibData();

frac = [1 .5 .1 .05];
scales = [.25 .5 1 2];
for i = 1:length(scales)
  opts.frac = frac(i);
  opts.scale = scales(i);
  [X,Y] = genCalibData(model,opts);
  model.beta(i) = learnBeta(X(:),Y(:));
end
