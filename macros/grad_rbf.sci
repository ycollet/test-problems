function [z]=grad_rbf(x)
  z=[ -2*(1-x(1))-420*x(1).*(x(2)-x(1).^2) ; 210*(x(2)-x(1).^2) ];
endfunction

