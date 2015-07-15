function W = LtoW( L )
% W = LtoW( L )

[h,w] = size(L);
W = zeros(h,w,121);

if(1)
  % simple version
  for y = 1:h
    for x = 1:w
      W(y,x,L(y,x)) = 1;
    end
  end
else
  % harder to read but faster
  [x,y] = meshgrid(1:w,1:h);
  x = x(:); y = y(:); z = L(:);
  W( y + (x-1)*h + (z-1)*h*w ) = 1;
end
