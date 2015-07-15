function Ipad = impad4( I, r )
% Pad image by r on all sides, plus extra on
% the right and bottom as necessary to make
% divisible by 4, as required by edgesChns
if nargin<2, r=0; end
[h,w,~] = size(I); p = [r r r r];
p([2 4]) = p([2 4]) + mod(4-mod([h w]+2*r,4),4);
Ipad = imPad( I, p, 'symmetric' );
end
