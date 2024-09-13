function Result = more_thuente_2(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("more_thuente_2: argument must be a scalar");
  end
  _beta  = 0.004;
  Result = (x + _beta)^5 - 2*(x + _beta)^4
endfunction
