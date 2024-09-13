ilib_for_link(['dmsafg','dmsasp','dmsahs'],['dmsafg.o','dmsasp.o','dmsahs.o'],[],'f');

// load the shared library 
exec loader.sce 

nx   = 4;
ny   = 4;
c    = 25.0;
bottom = ones(nx+2,1);
top    = ones(nx+2,1);
left   = ones(ny+2,1);
right  = ones(ny+2,1);

// Compute  function and gradient
// nx, ny: number of interior grid points
// c: angle of twist (typical value c=5)
//subroutine dmsafg(nx,ny,x,f,fgrad,task,bottom,top,left,right)
//character*(*) task
//integer nx, ny
//double precision f, bottom, top, left, right
//double precision x(nx*ny), fgrad(nx*ny)
// Compute Upper bound
task = 'XU';
x_upper = fort('dmsafg',nx,1,'i',ny,2,'i',task,6,'c',bottom,7,'d',top,8,'d',left,9,'d',right,10,'d','out',[nx*ny,1],3,'d');
// Compute Lower bound
task = 'XL';
x_lower = fort('dmsafg',nx,1,'i',ny,2,'i',task,6,'c',bottom,7,'d',top,8,'d',left,9,'d',right,10,'d','out',[nx*ny,1],3,'d');
// Compute Starting Point
task = 'XS';
x_start = fort('dmsafg',nx,1,'i',ny,2,'i',task,6,'c',bottom,7,'d',top,8,'d',left,9,'d',right,10,'d','out',[nx*ny,1],3,'d');
// Compute function and gradient
task = 'FG';
[f, fgrad] = fort('dmsafg',nx,1,'i',ny,2,'i',x_start,3,'d',task,6,'c',bottom,7,'d',top,8,'d',left,9,'d',right,10,'d','out',[1,1],4,'d',[nx*ny,1],5,'d');
printf('Result of dmsafg: f = %f\n',f);
disp(fgrad');
printf('Upper bound:'); disp(x_upper');
printf('Lower bound:'); disp(x_lower');
printf('Starting point:'); disp(x_start');

// Compute the sparsity structure of the hessian matrix
//subroutine dmsasp(nx,ny,nnz,indrow,indcol)
//integer nx, ny, nnz
//integer indrow(*), indcol(*)
indrow = [];
indcol = [];
[_nnz] = fort('dmsasp',nx,1,'i',ny,2,'i',indrow,4,'i',indcol,5,'i','out',[1,1],3,'i');
[indrow, indcol] = fort('dmsasp',nx,1,'i',ny,2,'i','out',[_nnz,1],4,'i',[_nnz,1],5,'i');
printf('Result of dmsasp: %d\n',_nnz);
disp(indrow'); disp(indcol');

// Compute y = H.s
//subroutine dmsahs(nx,ny,s,y,bottom,top,left,right)
//integer nx, ny
//double precision s(nx*ny), y(nx*ny)
// x -> x((j-1)*nx*i)
[y] = fort('dmsahs',nx,1,'i',ny,2,'i',x_start,3,'d',bottom,5,'d',top,6,'d',left,7,'d',right,8,'d','out',[nx*ny,1],4,'d');
printf('Result of dmsahs:'); 
disp(y');
