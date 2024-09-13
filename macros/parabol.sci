function [z]=parabol(x,y)
  [nargout, nargin] = argn();
  if (nargin==1) then
    z = x(1).^2 + x(2).^2;
  else
    z = x.^2 + y.^2;
  end
endfunction
