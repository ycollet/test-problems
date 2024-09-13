test_america          = %t;
test_bratu3db         = %f;
test_condor           = %f;
test_contor2          = %f;
test_contor           = %f;
test_ellipsoid        = %f;
test_genpack_cc_mina  = %f;
test_genpack_csq_mina = %f;
test_hardcube         = %f;
test_hardspheres      = %f;
test_kissing2         = %f;
test_kissing          = %f;
test_location         = %f;
test_mountain1        = %f;
test_mountain2        = %f
test_packccmn         = %f;
test_packccmn_feas    = %f;
test_packcrmn_feas    = %f;
test_pedazos4         = %f;
test_piecefit         = %f;
test_simfock2         = %f;
test_simfock          = %f;

//
// Test of the 'america' optimization problem
//

if test_america then
  // Test of the inip function
  printf('+------------------+\n');
  printf('| problem: america |\n');
  printf('+------------------+\n\n');
  
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('america');
  printf('inip: number of variables: %d - number of constraints: %d\n',n,m);
  printf('inip: vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('inip: vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('inip: initial point:'); disp(x(1:min(10,n)));
  printf('inip: vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('inip: vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('inip: vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('inip: vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
  // Test of the evalf function
  [f,flag] = evalf('america',n,x);
  printf('evalf: value of the objective function: %f - flag: %d\n',f,flag);
  // Test of the evalg function
  [g,flag] = evalg('america',n,x);
  printf('evalg: value of the gradient:'); disp(g);
  printf('evalg: flag: %d\n',flag);
  // Test of the evalh function
  [hnnz,hlin,hcol,hval,flag] = evalh('america',n,x);
  printf('evalh: number of non-empty values in the hessian matrix: %d\n',hnnz);
  printf('evalh: line coordinates:'); disp(hlin(1:min(10,hnnz)));
  printf('evalh: column coordinates:'); disp(hcol(1:min(10,hnnz)));
  printf('evalh: corresponding values:'); disp(hval(1:min(10,hnnz)));
  printf('evalh: flag:%d\n',flag);
  // Test of the evalc function
  [c,flag] = evalc('america',n,x,1);
  printf('evalc: value of the constraint %d: %f - flag: %d\n',1,c,flag);
  // Test of the evaljac function
  [jcnnz,jcvar,jcval,flag] = evaljac('america',n,x,1);
  printf('evaljac: number of non-empty values in the jacobian vector of constraint %d: %d\n',1,jcnnz);
  printf('evaljac: index coordinates:'); disp(jcvar(1:min(10,jcnnz)));
  printf('evaljac: corresponding values:'); disp(jcval(1:min(10,jcnnz)));
  printf('evaljac: flag:%d\n',flag);
  // Test of the evalhc function
  [hcnnz,hclin,hccol,hcval,flag] = evalhc('america',n,x,1);
  printf('evalhc: number of non-empty values in the hessian matrix: %d\n',hcnnz);
  printf('evalhc: line coordinates:'); disp(hclin(1:min(10,hcnnz)));
  printf('evalhc: column coordinates:'); disp(hccol(1:min(10,hcnnz)));
  printf('evalhc: corresponding values:'); disp(hcval(1:min(10,hcnnz)));
  printf('evalhc: flag:%d\n',flag);
  // Test of the evalhlp function
  [hp,goth,flag] = evalhlp('america',n,x,m,lambda,x,1);
  printf('evalhlp: hessian product vector:'); disp(hp(1:min(10,n)));
  printf('evalhlp: goth = %d - flag = %d\n',goth,flag);
  // Test of the endp function
  endp('america',n,x,l,u,m,lambda,rho,equatn,linear);
  printf('endp: success\n');
end

//
// Test of the 'bratu3db' optimization problem
//

if test_bratu3db then
  n = 10;
  seed = 153.7;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('bratu3db',n,seed);
  printf('problem: bratu3db\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'condor' optimization problem
//

if test_condor then
  n_in = 10;
  p_in = 100;
  q_in = 100;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('condor',n_in,p_in,q_in);
  printf('problem: condor\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'contor2' optimization problem
//

if test_contor2 then
  n_in = 10;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('contor2',n_in);
  printf('problem: contor2\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'contor' optimization problem
//

if test_contor then
  n_in = 10;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('contor',n_in);
  printf('problem: contor\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'ellipsoid' optimization problem
//

if test_ellipsoid then
  nd_in = 10;
  np_in = 100;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('ellispoid',nd_in,np_in);
  printf('problem: ellipsoid\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'genpack-cc-mina' optimization problem
//

if test_genpack_cc_mina then
  iterad_in = 0.01;
  nite_in   = 10;
  seed_in   = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('genpack-cc-mina',iterad_in,nite_in,seed_in);
  printf('problem: genpack_cc_mina\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'genpack-csq-mina' optimization problem
//

if test_genpack_csq_mina then
  iterad_in = 0.01;
  nite_in   = 10;
  seed_in   = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('genpack-csq-mina',iterad_in,nite_in,seed_in);
  printf('problem: genpack_csq_mina\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'hardcube' optimization problem
//

if test_hardcube then
  nd_in   = 10;
  np_in   = 100;
  seed_in = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('hardcube',nd_in,np_in,seed_in);
  printf('problem: hardcube\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'hardspheres' optimization problem
//

if test_hardspheres then
  nd_in   = 10;
  np_in   = 100;
  seed_in = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('hardspheres',nd_in,np_in,seed_in);
  printf('problem: hardspheres\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'kissing2' optimization problem
//

if test_kissing2 then
  nd_in   = 10;
  np_in   = 100;
  seed_in = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('kissing2',nd_in,np_in,seed_in);
  printf('problem: kissing2\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'kissing' optimization problem
//

if test_kissing then
  nd_in   = 10;
  np_in   = 100;
  seed_in = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('kissing',nd_in,np_in,seed_in);
  printf('problem: kissing\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'location' optimization problem
//

//if test_location then
//  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('location');
//  printf('problem: location\n');
//  printf('number of variables: %d - number of constraints: %d\n',n,m);
//  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
//  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
//  printf('initial point:'); disp(x(1:min(10,n)));
//  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
//  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
//  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
//  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
//end

//
// Test of the 'mountain1' optimization problem
//

if test_mountain1 then
  np_in    = 10;
  dmax2_in = 1.0;
  seed_in  = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('mountain1',np_in,dmax2_in,seed_in);
  printf('problem: mountain1\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'mountain2' optimization problem
//

if test_mountain2 then
  np_in    = 10;
  dmax2_in = 1.0;
  seed_in  = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('mountain1',np_in,dmax2_in,seed_in);
  printf('problem: mountain2\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'packccmn' optimization problem
//

if test_packccmn then
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('packccmn');
  printf('problem: packccmn\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'packccmn-feas' optimization problem
//

if test_packccmn_feas then
  nd_in     = 10;
  nite_in   = 100;
  iterad_in = 1.0;
  objrad_in = 1.0;
  seed_in   = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('packccmn-feas',nd_in,nite_in,iterad_in,objrad_in,seed_in);
  printf('problem: packccmn-feas\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'packcrmn-feas' optimization problem
//
if test_packcrmn_feas then
  nd_in     = 10;
  nite_in   = 100;
  iterad_in = 1.0;
  objdim_in = 1.0*ones(nd_in,1);
  seed_in   = 157.3;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('packcrmn-feas',nd_in,nite_in,iterad_in,objrad_in,seed_in);
  printf('problem: packcrmn-feas\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'pedazos4' optimization problem
//

if test_pedazos4 then
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('pedazos4');
  printf('problem: pedazos4\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'piecefit' optimization problem
//

if test_piecefit then
  n_in = 10;
  p_in = 20;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('piecefit',n_in,k_in);
  printf('problem: piecefit\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'simfock2' optimization problem
//

if test_simfock2 then
  n_in = 10;
  k_in = 20;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('simfock2',n_in,k_in);
  printf('problem: simfock2\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end

//
// Test of the 'simfock' optimization problem
//

if test_simfock then
  n_in = 10;
  k_in = 20;
  [n, x, l, u, m, lambda, rho, equatn, linear] = inip('simfock',n_in,k_in);
  printf('problem: simfock\n');
  printf('number of variables: %d - number of constraints: %d\n',n,m);
  printf('vector of upper bounds:'); disp(u(1:min(10,n)));
  printf('vector of lower bounds:'); disp(l(1:min(10,n)));
  printf('initial point:'); disp(x(1:min(10,n)));
  printf('vector of lagrange multipliers:'); disp(lambda(1:min(10,m)));
  printf('vector of initial penalty parameters:'); disp(rho(1:min(10,m)));
  printf('vector of equality constraints flag (1 = equality constraint):'); disp(equatn(1:min(10,m)));
  printf('vector of linear constraints flag (1 = linear constraint):'); disp(linear(1:min(10,m)));
end
