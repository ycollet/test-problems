function Result = grad_more_thuente_3(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("grad_more_thuente_3: argument must be a scalar");
  end
  _beta  = 0.01;
  l      = 39;
  if (x<=1-_beta) then
    phi0 = -1;
  elseif (x>=1+_beta) then
    phi0 = 1;
  else
    phi0 = 1.0/_beta*(x-1);
  end
  Result = phi0 + (1-_beta)*cos(l*%pi/2.0*x);
endfunction
