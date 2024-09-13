// CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
// nonlinear function minimization. To be used under the terms of the
// GNU General Public License (http://www.gnu.org/copyleft/gpl.html).
//
//
// Author: Nikolaus Hansen, 2003. 
// e-mail: hansen[at]bionik.tu-berlin.de
// URL: http://www.bionik.tu-berlin.de/user/niko
// References: See end of file. Last change: October, 27, 2004

// Modification for scilab :
// Author: Yann COLLETTE, 2006. 
// e-mail: yann[dot]colletet[at]renault[dot]com
// URL: http://ycollette.free.fr

////////////////////
// Test functions //
////////////////////

function f=fsphere(x)
f=sum(x.^2);
endfunction

function f=fschwefel(x)
f = 0;
for i = 1:size(x,1),
  f = f+sum(x(1:i))^2;
end
endfunction

function f=fcigar(x)
f = x(1)^2 + 1e6*sum(x(2:$).^2);
endfunction

function f=fcigtab(x)
f = x(1)^2 + 1e8*x($)^2 + 1e4*sum(x(2:($-1)).^2);
endfunction

function f=ftablet(x)
f = 1e6*x(1)^2 + sum(x(2:$).^2);
endfunction

function f=felli(x)
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
f=1e6.^((0:N-1)/(N-1)) * x.^2;
endfunction

function f=felli100(x)
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
f=1e4.^((0:N-1)/(N-1)) * x.^2;
endfunction

function f=fplane(x)
f=x(1);
endfunction

function f=ftwoaxes(x)
f = sum(x(1:floor($/2)).^2) + 1e6*sum(x(floor(1+$/2):$).^2);
endfunction

function f=fparabR(x)
f = -x(1) + 100*sum(x(2:$).^2);
endfunction

function f=fsharpR(x)
f = -x(1) + 100*norm(x(2:$));
endfunction

function f=frosen(x)
if (size(x,1) < 2) then error('dimension must be greater than one'); end
f = 100*sum((x(1:$-1).^2 - x(2:$)).^2) + sum((x(1:$-1)-1).^2);
endfunction

function f=fdiffpow(x)
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
f=sum(abs(x).^(2+10*(0:N-1)'/(N-1)));
endfunction

function f=frastrigin10(x)
N = size(x,1); if (N < 2) then error('dimension must be greater than one'); end
scale=10.^((0:N-1)'/(N-1));
f = 10*size(x,1) + sum((scale.*x).^2 - 10*cos(2*%pi*(scale.*x)));
endfunction

function f=frand(x)
f=rand(1,1);
endfunction

