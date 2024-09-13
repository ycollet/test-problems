function Result = grad_more_thuente_2(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("grad_more_thuente_2: argument must be a scalar");
  end
  _beta  = 0.004;
  Result = 5*(x + _beta)^4 - 8*(x + _beta)^3;
endfunction
