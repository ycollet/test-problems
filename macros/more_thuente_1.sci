function Result = more_thuente_1(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("more_thuente_1: argument must be a scalar");
  end
  _beta  = 2;
  Result = - x / (x^2 + _beta);
endfunction
