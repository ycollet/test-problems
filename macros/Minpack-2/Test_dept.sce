listoffiles = ['daerfj','datrfj','dchqfj','dcpffj','dcprfj','dctsfj','dedffj','deptfg','depths','deptsp','dfdcfj','dfdcjs','dfdcsp','dficfj','dficjs','dficsp', ...
               'dgdffj','dgl1fg','dgl1hs','dgl1sp','dgl2co','dgl2fg','dgl2hs','dgl2sp','dhhdfj','diacfj','diadfj','diaofj','diarfj','dierfj','dierjs','diersp', ...
               'dljcfg','dmsabc','dmsafg','dmsahs','dmsasp','dodcfg','dodchs','dodcps','dodcsp','dpjbds','dpjbfg','dpjbhs','dpjbsp','dsfdfj','dsfdjs','dsfdsp', ...
               'dsfifj','dsfijs','dsfisp','dsscfg','dsschs','dsscsp','saerfj','satrfj','schqfj','scpffj','scprfj','sctsfj','sedffj','septfg','sepths','septsp', ...
               'sfdcfj','sfdcjs','sfdcsp','sficfj','sficjs','sficsp','sgdffj','sgl1fg','sgl1hs','sgl1sp','sgl2co','sgl2fg','sgl2hs','sgl2sp','shhdfj','siacfj', ...
               'siadfj','siaofj','siarfj','sierfj','sierjs','siersp','sljcfg','smsabc','smsafg','smsahs','smsasp','sodcfg','sodchs','sodcps','sodcsp','spjbds', ...
               'spjbfg','spjbhs','spjbsp','ssfdfj','ssfdjs','ssfdsp','ssfifj','ssfijs','ssfisp','ssscfg','ssschs','ssscsp'];

if ~isdef('daerfj') then
  id = link('./libminpack-2.so',listoffiles);
end

nx   = 4;
ny   = 4;
c    = 5.0;

// Compute  function and gradient
// nx, ny: number of interior grid points
// c: angle of twist (typical value c=5)
//subroutine deptfg(nx,ny,x,f,fgrad,task,c)
//character*(*) task
//integer nx, ny
//double precision f, c
//double precision x(nx*ny), fgrad(nx*ny)
// Compute Upper bound
task = 'XU';
x_upper = call('deptfg',nx,1,'i',ny,2,'i',task,6,'c',c,7,'d','out',[nx*ny,1],3,'d');
// Compute Lower bound
task = 'XL';
x_lower = call('deptfg',nx,1,'i',ny,2,'i',task,6,'c',c,7,'d','out',[nx*ny,1],3,'d');
// Compute Starting Point
task = 'XS';
x_start = call('deptfg',nx,1,'i',ny,2,'i',task,6,'c',c,7,'d','out',[nx*ny,1],3,'d');
// Compute function and gradient
task = 'FG';
[f, fgrad] = call('deptfg',nx,1,'i',ny,2,'i',x_start,3,'d',task,6,'c',c,7,'d','out',[1,1],4,'d',[nx*ny,1],5,'d');
printf('Result of deptfg: f = %f\n',f);
disp(fgrad');
printf('Upper bound:'); disp(x_upper');
printf('Lower bound:'); disp(x_lower');
printf('Starting point:'); disp(x_start');

// Compute the sparsity structure of the hessian matrix
//subroutine deptsp(nx,ny,nnz,indrow,indcol)
//integer nx, ny, nnz
//integer indrow(*), indcol(*)
indrow = [];
indcol = [];
[_nnz] = call('deptsp',nx,1,'i',ny,2,'i',indrow,4,'i',indcol,5,'i','out',[1,1],3,'i');
[indrow, indcol] = call('deptsp',nx,1,'i',ny,2,'i','out',[_nnz,1],4,'i',[_nnz,1],5,'i');
printf('Result of deptsp: %d\n',_nnz);
disp(indrow'); disp(indcol');

// Compute y = H.s
//subroutine depths(nx,ny,s,y)
//integer nx, ny
//double precision s(nx*ny), y(nx*ny)
// x -> x((j-1)*nx*i)
[y] = call('depths',nx,1,'i',ny,2,'i',x_start,3,'d','out',[nx*ny,1],4,'d');
printf('Result of depths:'); 
//disp(y');
