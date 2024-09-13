function Result = more_thuente_3(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("more_thuente_3: argument must be a scalar");
  end
  _beta  = 0.01;
  l      = 39;
  if (x<=1-_beta) then
    phi0 = 1 - x;
  elseif (x>=1+_beta) then
    phi0 = x - 1;
  else
    phi0 = 1.0/(2.0*_beta)*(x-1)^2+ 1.0/2.0 * _beta;
  end
  Result = phi0 + 2*(1-_beta)/(l*%pi)*sin(l*%pi/2.0*x);
endfunction
