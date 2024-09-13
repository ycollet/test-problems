function plotlandscape(upper,lower,N)

//--------------------------------------------------------------------------
//This is the plotting function of the Gaussian landscape generator(2D only)
//
//Syntax: plotlandscape(upper,lower,N)
//
//Example:plotlandscape(5,-5,100)
//
//Inputs:
//       upper boundary of the search space
//       lower boundary of the search space
//       number of sampling points
//
//Outputs:
//       3D surface plot
//       2D contour plot
//
//Author: Bo Yuan (boyuan@itee.uq.edu.au)
//--------------------------------------------------------------------------

[nargout,nargin] = argn();

if nargin~=3 then
  error('plotlandscape: Usage: plotlandscape(upper,lower,N)');
end

if upper<=lower then
  error('plotlandscape: Upper must be larger than Lower!');
end

if N<10 then
  error('plotlandscape: Incorrect N value!');
end

//-----------------------------------------------------

inc=(upper-lower)/N;

x=lower:inc:upper;     //x coordinates
y=lower:inc:upper;     //y coordinates

for i=1:length(x)
  for j=1:length(y)
    X(i,j) = x(i);
    Y(i,j) = y(j);
    Z(i,j) = fitness([x(i) y(j)]);
  end
end

scf;
drawlater;
subplot(2,1,1);
Color = graycolormap(128); //3D surface plot
surf(X,Y,Z);
xset('colormap',Color(64:$,:));
xtitle('Gaussian landscape','x1','x2','x3');
subplot(2,1,2);
contour(x,y,Z,[0.1 0.3 0.6 0.9]);  //2D contour plot
xtitle('Contour of the gaussian landscape','x1','x2');
drawnow;
endfunction
