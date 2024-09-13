function Result = YOK_3(x)
  [nargout, nargin] = argn();
  if (nargin>1) then
    error("YOK_3: argument must be a scalar");
  end
  _beta1 = 0.001;
  _beta2 = 0.01;
  Result = YOK(x,_beta1,_beta2);
endfunction
