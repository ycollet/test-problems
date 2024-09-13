function Result = grad_more_thuente_1(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("grad_more_thuente_1: argument must be a scalar");
  end
  _beta  = 2;
  Result = (x^2 - _beta) / (x^2 + _beta)^2;
endfunction
