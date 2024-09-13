//
//--------------------Lennard-Johns functions (lj and ljgrad computation)-----------------------------------------
//
function f = lj(x)                         // lj function
n = max (size (x))/ 3;
f = 0;
for i = 1 : n
   for j = i+1 : n
      i1 = 3*(i-1) + 1; i2 = i1 + 1; i3 = i2 + 1;
      j1 = 3*(j-1) + 1; j2 = j1 + 1; j3 = j2 + 1;
      r1 = x(i1)-x(j1); r2 = x(i2)-x(j2); r3 = x(i3)-x(j3);
      rr = r1 * r1 + r2 * r2 + r3 * r3;
      r6 = 1 / rr / rr / rr;
      f = f + (r6 - 2) * r6;
   end
end
endfunction
//
function g = ljgrad(x)                     // gradient of lj function
n = max (size (x)) / 3;
g = zeros(1,3*n);
for i = 1 : n
   for j = i+1 : n
      i1 = 3*(i-1) + 1; i2 = i1 + 1; i3 = i2 + 1;
      j1 = 3*(j-1) + 1; j2 = j1 + 1; j3 = j2 + 1;
      r1 = x(i1)-x(j1); r2 = x(i2)-x(j2); r3 = x(i3)-x(j3);
      rr = r1 * r1 + r2 * r2 + r3 * r3;
      r = sqrt (rr);
      r6 = 1 / rr / rr / rr;
      dr = - 12 * (r6 - 1) * r6 / r;
      g(i1) = g(i1) + dr * r1 / r;
      g(i2) = g(i2) + dr * r2 / r;
      g(i3) = g(i3) + dr * r3 / r;
      g(j1) = g(j1) - dr * r1 / r;
      g(j2) = g(j2) - dr * r2 / r;
      g(j3) = g(j3) - dr * r3 / r;
   end
end
endfunction

function plot_lj(x)
  for j=1:(N-1)
    for k=j+1:N
      plot3d([x(3*j-2);x(3*k-2)],[x(3*j-1);x(3*k-1)],[x(3*j);x(3*k)]);
      plot3d([x(3*j-2);x(3*k-2)],[x(3*j-1);x(3*k-1)],[x(3*j);x(3*k)],-1);
    end
  end
endfunction


























 
     

