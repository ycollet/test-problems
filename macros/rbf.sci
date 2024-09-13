function [z]=rbf(x,y)
  [nargout, nargin] = argn();
  if (nargin==1) then
    z = (1-x(1))^2+105*(x(2)-x(1)^2)^2;
  else
    z = (1-x).^2+105*(y-x.^2).^2;
  end
endfunction

