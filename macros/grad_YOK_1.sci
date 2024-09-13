function Result = grad_YOK_1(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("grad_YOK_1: argument must be a scalar");
  end
  _beta1 = 0.001;
  _beta2 = 0.001;
  Result = grad_YOK(x, _beta1, _beta2);
endfunction
