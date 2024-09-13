function Result = grad_YOK(x,_beta1,_beta2)
  [nargout, nargin] = argn();
  if ((size(x,1)>1)|(size(x,2)>1)) then
    error("grad_YOK: argument must be a scalar");
  end
  _lambda1 = (1+_beta1^2)^0.5 - _beta1;
  _lambda2 = (1+_beta2^2)^0.5 - _beta2;
  Result = - 0.5*_lambda1*2*(1-x)*((1-x)^2+_beta2^2)^(-0.5) - 0.5*_lambda2*2*x*(x^2+_beta1^2)^(-0.5);
endfunction
