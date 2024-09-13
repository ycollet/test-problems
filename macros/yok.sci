function Result = YOK(x, _beta1, _beta2)
  [nargout, nargin] = argn();
  if ((size(x,1)>1)|(size(x,2)>1)) then
    error("YOK: argument must be a scalar");
  end
  _lambda1 = (1+_beta1^2)^0.5 - _beta1;
  _lambda2 = (1+_beta2^2)^0.5 - _beta2;
  Result = _lambda1*((1-x)^2+_beta2^2)^0.5 + _lambda2*(x^2+_beta1^2)^0.5;
endfunction
