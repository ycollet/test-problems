function Result = YOK_2(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("YOK_2: argument must be a scalar");
  end
  _beta1 = 0.01;
  _beta2 = 0.001;
  Result = YOK(x,_beta1,_beta2);
endfunction
