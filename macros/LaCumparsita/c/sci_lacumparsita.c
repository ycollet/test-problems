#include <stack-c.h>
#include <string.h>
#include <stdio.h>

#define dstk(x) (((double *)C2F(stack).Stk) + x - 1)
#define lstk(x) (((logical *)C2F(stack).Stk) + x - 1)

#include "lacumparsita.h"

////////////////////
// inip interface //
////////////////////

int interface_inip_(char * fname)
{
  static int minlhs = 9, maxlhs = 9, minrhs = 1, maxrhs = 3;
  int l1 = 0, m1, n1; // variable n
  int l2 = 0, m2, n2; // variable x(n)
  int l3 = 0, m3, n3; // variable l(n)
  int l4 = 0, m4, n4; // variable u(n)
  int l5 = 0, m5, n5; // variable m
  int l6 = 0, m6, n6; // variable lambda(m)
  int l7 = 0, m7, n7; // variable rho(m)
  int l8 = 0, m8, n8; // variable equatn(m)
  int l9 = 0, m9, n9; // variable linear(m)
  int l10 = 0, m10, n10; // variable function name
  int l11 = 0, m11, n11; // input variable 1
  int l12 = 0, m12, n12; // input variable 2
  int l13 = 0, m13, n13; // input variable 3
  int l14 = 0, m14, n14; // input variable 4
  int l15 = 0, m15, n15; // input variable 5
  int j,i, Index = 1;
  int ierr = 0;

  integer    * out_n       = NULL;
  doublereal * out_xn      = NULL;
  doublereal * out_ln      = NULL;
  doublereal * out_un      = NULL;
  integer    * out_m       = NULL;
  doublereal * out_lambdam = NULL;
  doublereal * out_rhom    = NULL;
  logical    * out_equatnm = NULL;
  logical    * out_linearm = NULL;

  GetRhsVar(1,"c",&m10, &n10, &l10); // Variable function name

  sciprint("name = **%s**\n",cstk(l10));

  if (strcmp(cstk(l10),"america")==0)
    {
      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(264*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(264*sizeof(doublereal));
      out_un      = (doublereal *)malloc(264*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(17*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(17*sizeof(doublereal));
      out_equatnm = (logical *)malloc(17*sizeof(logical));
      out_linearm = (logical *)malloc(17*sizeof(logical));

      ierr = america_inip__(out_n,        // Variable n         - output
			    out_xn,       // Variable x(n)      - output
			    out_ln,       // Variable l(n)      - output
			    out_un,       // Variable u(n)      - output
			    out_m,        // Variable m         - output
			    out_lambdam,  // Variable lambda(m) - output
			    out_rhom,     // Variable rho(m)    - output
			    out_equatnm,  // Variable equatn(m) - output
			    out_linearm); // Variable linear(m) - output
    }
  else if (strcmp(cstk(l10),"bratu3db")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable np
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable seed

      Index = 3;

      int    n    = *istk(l11); n = (n>100)?100:n;
      double seed = *istk(l12);

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n^3*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n^3*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n^3*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc((n-2)^3*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc((n-2)^3*sizeof(doublereal));
      out_equatnm = (logical *)malloc((n-2)^3*sizeof(logical));
      out_linearm = (logical *)malloc((n-2)^3*sizeof(logical));

      ierr = bratu3db_inip__(out_n,        // Variable n         - output
			     out_xn,       // Variable x(n)      - output
			     out_ln,       // Variable l(n)      - output
			     out_un,       // Variable u(n)      - output
			     out_m,        // Variable m         - output
		 	     out_lambdam,  // Variable lambda(m) - output
			     out_rhom,     // Variable rho(m)    - output
			     out_equatnm,  // Variable equatn(m) - output
			     out_linearm,  // Variable linear(m) - output
			     &n,           // Variable np        - input
			     &seed);       // Variable seed      - input
    }
  else if (strcmp(cstk(l10),"cache")==0)
    {
      int n = 3*128;
      int m = 2*128+2;

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = cache_inip__(out_n,        // Variable n         - output
			  out_xn,       // Variable x(n)      - output
			  out_ln,       // Variable l(n)      - output
			  out_un,       // Variable u(n)      - output
			  out_m,        // Variable m         - output
		 	  out_lambdam,  // Variable lambda(m) - output
			  out_rhom,     // Variable rho(m)    - output
			  out_equatnm,  // Variable equatn(m) - output
			  out_linearm); // Variable linear(m) - output

    }
  else if (strcmp(cstk(l10),"condor")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable n_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable p_in
      GetRhsVar(4,"i",&m13, &n13, &l13); // Variable q_in

      Index = 4;

      int n = *istk(l11); // n = n_in
      int m = *istk(l11); // m = n - 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = condor_inip__(out_n,        // Variable n         - output
			   out_xn,       // Variable x(n)      - output
			   out_ln,       // Variable l(n)      - output
			   out_un,       // Variable u(n)      - output
			   out_m,        // Variable m         - output
		 	   out_lambdam,  // Variable lambda(m) - output
			   out_rhom,     // Variable rho(m)    - output
			   out_equatnm,  // Variable equatn(m) - output
			   out_linearm,  // Variable linear(m) - output
			   istk(l11),    // Variable n_in      - input
			   istk(l12),    // Variable p_in      - input
			   istk(l13));   // Variable q_in      - input
    }
  else if (strcmp(cstk(l10),"contor2")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable n_in

      Index = 2;

      int n = *istk(l11) + 3; // n = n_in + 3
      int m = *istk(l11) - 2; // m = n - 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = contor2_inip__(out_n,        // Variable n         - output
			    out_xn,       // Variable x(n)      - output
			    out_ln,       // Variable l(n)      - output
			    out_un,       // Variable u(n)      - output
			    out_m,        // Variable m         - output
		 	    out_lambdam,  // Variable lambda(m) - output
			    out_rhom,     // Variable rho(m)    - output
			    out_equatnm,  // Variable equatn(m) - output
			    out_linearm,  // Variable linear(m) - output
			    istk(l11));   // Variable n_in      - input
    }
  else if (strcmp(cstk(l10),"contor")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable n_in

      Index = 2;

      int n = *istk(l11); // n = n_in
      int m = *istk(l11); // m = n - 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = contor_inip__(out_n,        // Variable n         - output
			   out_xn,       // Variable x(n)      - output
			   out_ln,       // Variable l(n)      - output
			   out_un,       // Variable u(n)      - output
			   out_m,        // Variable m         - output
		 	   out_lambdam,  // Variable lambda(m) - output
			   out_rhom,     // Variable rho(m)    - output
			   out_equatnm,  // Variable equatn(m) - output
			   out_linearm,  // Variable linear(m) - output
			   istk(l11));   // Variable n_in      - input
    }
  else if (strcmp(cstk(l10),"ellipsoid")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable nd_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable np_in

      Index = 3;

      int n = *istk(l11) * (*istk(l11) + 1) / 2; // n = nd_in*(nd_in-1)/2
      int m = *istk(l12); // m = np_in

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = ellipsoid_inip__(out_n,        // Variable n         - output
			      out_xn,       // Variable x(n)      - output
			      out_ln,       // Variable l(n)      - output
			      out_un,       // Variable u(n)      - output
			      out_m,        // Variable m         - output
		 	      out_lambdam,  // Variable lambda(m) - output
			      out_rhom,     // Variable rho(m)    - output
			      out_equatnm,  // Variable equatn(m) - output
			      out_linearm,  // Variable linear(m) - output
			      istk(l11),    // Variable nd_in     - input
			      istk(l12));   // Variable np_in     - input
    }
  else if (strcmp(cstk(l10),"genpack-cc-mina")==0)
    {
      GetRhsVar(2,"d",&m11, &n11, &l11); // Variable iterad_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable nite_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable seed_in

      Index = 4;

      int n = *istk(l12) * 2 + 1; // n = nite_in * 2 + 1 
      int m = *istk(l12) + 1; // m = nite_in + 1

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = genpack_cc_mina_inip__(out_n,        // Variable n         - output
			            out_xn,       // Variable x(n)      - output
			            out_ln,       // Variable l(n)      - output
			            out_un,       // Variable u(n)      - output
			            out_m,        // Variable m         - output
		 	            out_lambdam,  // Variable lambda(m) - output
			            out_rhom,     // Variable rho(m)    - output
			            out_equatnm,  // Variable equatn(m) - output
			            out_linearm,  // Variable linear(m) - output
			            dstk(l11),    // Variable iterad_in - input
			            istk(l12),    // Variable nite_in   - input
			            dstk(l13));   // Variable seed_in   - input
    }
  else if (strcmp(cstk(l10),"genpack-csq-mina")==0)
    {
      GetRhsVar(2,"d",&m11, &n11, &l11); // Variable iterad_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable nite_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable seed_in

      Index = 4;

      int n = *istk(l12) * 2 + 1; // n = nite_in * 2 + 1
      int m = *istk(l12) * 2 + 1; // m = nite_in * 2 + 1

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = genpack_csq_mina_inip__(out_n,        // Variable n         - output
			             out_xn,       // Variable x(n)      - output
			             out_ln,       // Variable l(n)      - output
			             out_un,       // Variable u(n)      - output
			             out_m,        // Variable m         - output
		 	             out_lambdam,  // Variable lambda(m) - output
			             out_rhom,     // Variable rho(m)    - output
			             out_equatnm,  // Variable equatn(m) - output
			             out_linearm,  // Variable linear(m) - output
			             dstk(l11),    // Variable iterad_in - input
			             istk(l12),    // Variable nite_in   - input
			             dstk(l13));   // Variable seed_in   - input
    }
  else if (strcmp(cstk(l10),"hardcube")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable nd_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable np_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable seed_in

      Index = 4;

      int n = *istk(l11) * *istk(l12) + 1; // n = nd_in * np_in + 1
      int m = *istk(l12) + *istk(l12) * (*istk(l12) - 1) / 2; // m = np_in + np_in*(np_in - 1) / 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = hardcube_inip__(out_n,        // Variable n         - output
			     out_xn,       // Variable x(n)      - output
			     out_ln,       // Variable l(n)      - output
			     out_un,       // Variable u(n)      - output
			     out_m,        // Variable m         - output
		 	     out_lambdam,  // Variable lambda(m) - output
			     out_rhom,     // Variable rho(m)    - output
			     out_equatnm,  // Variable equatn(m) - output
			     out_linearm,  // Variable linear(m) - output
			     istk(l11),    // Variable nd_in      - input
			     istk(l12),    // Variable np_in      - input
			     dstk(l13));   // Variable seed_in      - input
    }
  else if (strcmp(cstk(l10),"hardspheres")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable nd_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable np_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable seed_in

      Index = 4;

      int n = *istk(l11) * *istk(l12) + 1; // n = nd_in * nd_in + 1
      int m = *istk(l12) + *istk(l12) * (*istk(l12) - 1) / 2; // m = np_in + np_in * (np_in - 1) / 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = hardspheres_inip__(out_n,        // Variable n         - output
			        out_xn,       // Variable x(n)      - output
			        out_ln,       // Variable l(n)      - output
			        out_un,       // Variable u(n)      - output
			        out_m,        // Variable m         - output
		 	        out_lambdam,  // Variable lambda(m) - output
			        out_rhom,     // Variable rho(m)    - output
			        out_equatnm,  // Variable equatn(m) - output
			        out_linearm,  // Variable linear(m) - output
			        istk(l11),    // Variable nd_in      - input
			        istk(l12),    // Variable np_in      - input
			        dstk(l13));   // Variable seed_in      - input
    }
  else if (strcmp(cstk(l10),"kissing2")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable nd_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable np_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable seed_in

      Index = 4;

      int n = *istk(l11) * *istk(l12); // n = nd_in * np_in
      int m = *istk(l12) + *istk(l12) * (*istk(l12) - 1) / 2; // m = np_in + np_in  * (np_in -1)/2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = kissing2_inip__(out_n,        // Variable n         - output
			     out_xn,       // Variable x(n)      - output
			     out_ln,       // Variable l(n)      - output
			     out_un,       // Variable u(n)      - output
			     out_m,        // Variable m         - output
		 	     out_lambdam,  // Variable lambda(m) - output
			     out_rhom,     // Variable rho(m)    - output
			     out_equatnm,  // Variable equatn(m) - output
			     out_linearm,  // Variable linear(m) - output
			     istk(l11),    // Variable nd_in      - input
			     istk(l12),    // Variable np_in      - input
			     dstk(l13));   // Variable seed_in    - input
    }
  else if (strcmp(cstk(l10),"kissing")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable nd_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable np_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable seed_in

      Index = 4;

      int n = *istk(l11) * *istk(l12); // n = nd_in * np_in
      int m = *istk(l12) + *istk(l12) * (*istk(l12) - 1) / 2; // m = n - 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = kissing_inip__(out_n,        // Variable n         - output
			    out_xn,       // Variable x(n)      - output
			    out_ln,       // Variable l(n)      - output
			    out_un,       // Variable u(n)      - output
			    out_m,        // Variable m         - output
		 	    out_lambdam,  // Variable lambda(m) - output
			    out_rhom,     // Variable rho(m)    - output
			    out_equatnm,  // Variable equatn(m) - output
			    out_linearm,  // Variable linear(m) - output
			    istk(l11),    // Variable nd_in     - input
			    istk(l12),    // Variable np_in     - input
			    dstk(l13));   // Variable seed_in   - input
    }
  else if (strcmp(cstk(l10),"location")==0)
    {
      //
      // YC: a revoir: il y a des paramètres ...
      //
      int n = *istk(l11); // n = n_in
      int m = *istk(l11); // m = n - 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = location_inip__(out_n,        // Variable n         - output
			     out_xn,       // Variable x(n)      - output
			     out_ln,       // Variable l(n)      - output
			     out_un,       // Variable u(n)      - output
			     out_m,        // Variable m         - output
		 	     out_lambdam,  // Variable lambda(m) - output
			     out_rhom,     // Variable rho(m)    - output
			     out_equatnm,  // Variable equatn(m) - output
			     out_linearm); // Variable linear(m) - output
    }
  else if (strcmp(cstk(l10),"mountain1")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable np_in
      GetRhsVar(3,"d",&m12, &n12, &l12); // Variable dmax2_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable seed_in

      Index = 4;

      int n = 2 * *istk(l11) + 1; // n = 2 * np_in + 1
      int m = 2 * *istk(l11) + 1; // m = 2 * np_in + 1

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = mountain1_inip__(out_n,        // Variable n         - output
			      out_xn,       // Variable x(n)      - output
			      out_ln,       // Variable l(n)      - output
			      out_un,       // Variable u(n)      - output
			      out_m,        // Variable m         - output
		 	      out_lambdam,  // Variable lambda(m) - output
			      out_rhom,     // Variable rho(m)    - output
			      out_equatnm,  // Variable equatn(m) - output
			      out_linearm,  // Variable linear(m) - output
			      istk(l11),    // Variable np_in     - input
			      dstk(l12),    // Variable dmax2_in  - input
			      dstk(l13));   // Variable seed_in   - input
    }
  else if (strcmp(cstk(l10),"mountain2")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable np_in
      GetRhsVar(3,"d",&m12, &n12, &l12); // Variable dmax2_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable seed_in

      Index = 4;

      int n = 2 * *istk(l11) + 1; // n = 2 * np_in + 1
      int m = 2 * *istk(l11) + 1; // m = 2 * np_in + 1

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = mountain2_inip__(out_n,        // Variable n         - output
			      out_xn,       // Variable x(n)      - output
			      out_ln,       // Variable l(n)      - output
			      out_un,       // Variable u(n)      - output
			      out_m,        // Variable m         - output
		 	      out_lambdam,  // Variable lambda(m) - output
			      out_rhom,     // Variable rho(m)    - output
			      out_equatnm,  // Variable equatn(m) - output
			      out_linearm,  // Variable linear(m) - output
			      istk(l11),    // Variable np_in     - input
			      dstk(l12),    // Variable dmax2_in  - input
			      dstk(l13));   // Variable seed_in   - input
    }
  else if (strcmp(cstk(l10),"packccmn")==0)
    {
      int n = 2 * 80;
      int m = 80;    

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = packccmn_inip__(out_n,        // Variable n         - output
			     out_xn,       // Variable x(n)      - output
			     out_ln,       // Variable l(n)      - output
			     out_un,       // Variable u(n)      - output
			     out_m,        // Variable m         - output
		 	     out_lambdam,  // Variable lambda(m) - output
			     out_rhom,     // Variable rho(m)    - output
			     out_equatnm,  // Variable equatn(m) - output
			     out_linearm); // Variable linear(m) - output
    }
  else if (strcmp(cstk(l10),"packccmn-feas")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable nd_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable nite_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable iterad_in
      GetRhsVar(5,"d",&m14, &n14, &l14); // Variable objrad_in
      GetRhsVar(6,"d",&m15, &n15, &l15); // Variable seed_in

      Index = 6;

      int n = *istk(l11) * *istk(l12); // n = nd_in * nite_in
      int m = *istk(l12) + *istk(l12) * (*istk(l12) - 1) / 2; // m = nite_in + nite_in * (nite_in - 1) / 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = packccmn_feas_inip__(out_n,        // Variable n         - output
			          out_xn,       // Variable x(n)      - output
			          out_ln,       // Variable l(n)      - output
			          out_un,       // Variable u(n)      - output
			          out_m,        // Variable m         - output
		                  out_lambdam,  // Variable lambda(m) - output
			          out_rhom,     // Variable rho(m)    - output
			          out_equatnm,  // Variable equatn(m) - output
			          out_linearm,  // Variable linear(m) - output
			          istk(l11),    // Variable nd_in     - input
			          istk(l12),    // Variable nite_in   - input
			          dstk(l13),    // Variable iterad_in - input
			          dstk(l14),    // Variable objrad_in - input
			          dstk(l15));   // Variable seed_in   - input
    }
  else if (strcmp(cstk(l10),"packcrmn-feas")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable nd_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable nite_in
      GetRhsVar(4,"d",&m13, &n13, &l13); // Variable iterad_in
      GetRhsVar(5,"d",&m14, &n14, &l14); // Variable objdim_in - YC: verifier la concordance de la dimension pour ce parametre
      GetRhsVar(6,"d",&m15, &n15, &l15); // Variable seed_in

      Index = 6;

      int n = *istk(l11) * *istk(l12); // n = nd_in * nite_in
      int m = *istk(l12) * (*istk(l12) - 1) / 2; // m = nite_in * (nite_in - 1) / 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = packccmn_feas_inip__(out_n,        // Variable n         - output
			          out_xn,       // Variable x(n)      - output
			          out_ln,       // Variable l(n)      - output
			          out_un,       // Variable u(n)      - output
			          out_m,        // Variable m         - output
		                  out_lambdam,  // Variable lambda(m) - output
			          out_rhom,     // Variable rho(m)    - output
			          out_equatnm,  // Variable equatn(m) - output
			          out_linearm,  // Variable linear(m) - output
			          istk(l11),    // Variable nd_in     - input
			          istk(l12),    // Variable nite_in   - input
			          dstk(l13),    // Variable iterad_in - input
			          dstk(l14),    // Variable objdim_in - input
			          dstk(l15));   // Variable seed_in   - input
    }
  else if (strcmp(cstk(l10),"pedazos4")==0)
    {
      int n = 9;
      int m = 2;

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = pedazos4_inip__(out_n,        // Variable n         - output
			     out_xn,       // Variable x(n)      - output
			     out_ln,       // Variable l(n)      - output
			     out_un,       // Variable u(n)      - output
			     out_m,        // Variable m         - output
		             out_lambdam,  // Variable lambda(m) - output
			     out_rhom,     // Variable rho(m)    - output
			     out_equatnm,  // Variable equatn(m) - output
			     out_linearm); // Variable linear(m) - output
    }
  else if (strcmp(cstk(l10),"piecefit")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable n_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable p_in

      Index = 3;

      int n = *istk(l11); // n = n_in
      int m = 1;

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = piecefit_inip__(out_n,        // Variable n         - output
			     out_xn,       // Variable x(n)      - output
			     out_ln,       // Variable l(n)      - output
			     out_un,       // Variable u(n)      - output
			     out_m,        // Variable m         - output
		             out_lambdam,  // Variable lambda(m) - output
			     out_rhom,     // Variable rho(m)    - output
			     out_equatnm,  // Variable equatn(m) - output
			     out_linearm,  // Variable linear(m) - output
			     istk(l11),    // Variable n_in      - input
			     istk(l12));   // Variable p_in      - input
    }
  else if (strcmp(cstk(l10),"simfock2")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable n_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable k_in

      Index = 3;

      int n = *istk(l11) * *istk(l12); // n = n_in * k_in
      int m = *istk(l11) + *istk(l11) * (*istk(l11) - 1) / 2; // m = n_in + n_in * (n_in - 1) / 2

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = simfock2_inip__(out_n,        // Variable n         - output
			     out_xn,       // Variable x(n)      - output
			     out_ln,       // Variable l(n)      - output
			     out_un,       // Variable u(n)      - output
			     out_m,        // Variable m         - output
		             out_lambdam,  // Variable lambda(m) - output
			     out_rhom,     // Variable rho(m)    - output
			     out_equatnm,  // Variable equatn(m) - output
			     out_linearm,  // Variable linear(m) - output
			     istk(l11),    // Variable n_in      - input
			     istk(l12));   // Variable k_in      - input
    }
  else if (strcmp(cstk(l10),"simfock")==0)
    {
      GetRhsVar(2,"i",&m11, &n11, &l11); // Variable n_in
      GetRhsVar(3,"i",&m12, &n12, &l12); // Variable k_in

      Index = 3;

      int n = *istk(l12) * *istk(l12); // n = k_in * k_in
      int m = 2 * *istk(l12) * *istk(l12) + 1; // m = 2 * k_in * k_in + 1

      out_n       = (integer *)malloc(1*sizeof(integer));
      out_xn      = (doublereal *)malloc(n*sizeof(doublereal));
      out_ln      = (doublereal *)malloc(n*sizeof(doublereal));
      out_un      = (doublereal *)malloc(n*sizeof(doublereal));
      out_m       = (integer *)malloc(1*sizeof(integer));
      out_lambdam = (doublereal *)malloc(m*sizeof(doublereal));
      out_rhom    = (doublereal *)malloc(m*sizeof(doublereal));
      out_equatnm = (logical *)malloc(m*sizeof(logical));
      out_linearm = (logical *)malloc(m*sizeof(logical));

      ierr = simfock_inip__(out_n,        // Variable n         - output
			    out_xn,       // Variable x(n)      - output
			    out_ln,       // Variable l(n)      - output
			    out_un,       // Variable u(n)      - output
			    out_m,        // Variable m         - output
		            out_lambdam,  // Variable lambda(m) - output
			    out_rhom,     // Variable rho(m)    - output
			    out_equatnm,  // Variable equatn(m) - output
			    out_linearm,  // Variable linear(m) - output
			    istk(l11),    // Variable n_in      - input
			    istk(l12));   // Variable k_in      - input
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  // Variable n
  m1 = 1; n1 = 1;
  // Variable x(n)
  m2 = 1; n2 = *out_n;
  // Variable l(n)
  m3 = 1; n3 = *out_n;
  // Variable u(n)
  m4 = 1; n4 = *out_n;
  // Variable m
  m5 = 1; n5 = 1;
  // Variable lambda(m)
  m6 = 1; n6 = *out_m;
  // Variable rho(m)
  m7 = 1; n7 = *out_m;
  // Variable equatn(m)
  m8 = 1; n8 = *out_m;
  // Variable linear(m)
  m9 = 1; n9 = *out_m;

  CreateVar(Index + 1, "i",&m1, &n1, &l1); // Variable n
  CreateVar(Index + 2, "d",&m2, &n2, &l2); // Variable x(n)
  CreateVar(Index + 3, "d",&m3, &n3, &l3); // Variable l(n)
  CreateVar(Index + 4, "d",&m4, &n4, &l4); // Variable u(n)
  CreateVar(Index + 5, "i",&m5, &n5, &l5); // Variable m
  CreateVar(Index + 6, "d",&m6, &n6, &l6); // Variable lambda(m)
  CreateVar(Index + 7, "d",&m7, &n7, &l7); // Variable rho(m)
  CreateVar(Index + 8, "d",&m8, &n8, &l8); // Variable equatn(m)
  CreateVar(Index + 9, "d",&m9, &n9, &l9); // Variable linear(m)

  *istk(l1) = *out_n;
  *istk(l5) = *out_m;

  for(i=0;i<*out_n;i++)
    {
      *(dstk(l2)+i) = out_xn[i];
      *(dstk(l3)+i) = out_ln[i];
      *(dstk(l4)+i) = out_un[i];
    } /* End For */

  for(i=0;i<*out_m; i++)
    {
      *(dstk(l6)+i) = out_lambdam[i];
      *(dstk(l7)+i) = out_rhom[i];
      *(dstk(l8)+i) = out_equatnm[i];
      *(dstk(l9)+i) = out_linearm[i];
    } /* End For */
  
  free(out_n);
  free(out_xn);
  free(out_ln);
  free(out_un);
  free(out_m);
  free(out_lambdam);
  free(out_rhom);
  free(out_equatnm);
  free(out_linearm);

  LhsVar(1) = 2;
  LhsVar(2) = 3;
  LhsVar(3) = 4;
  LhsVar(4) = 5;
  LhsVar(5) = 6;
  LhsVar(6) = 7;
  LhsVar(7) = 8;
  LhsVar(8) = 9;
  LhsVar(9) = 10;

  PutLhsVar();

  return 0;
}

/////////////////////
// interface evalf //
/////////////////////


int interface_evalf_(char * fname)
{
  int l1, m1, n1; // variable n
  int l2, m2, n2; // variable x(n)
  int l3, m3, n3; // variable f
  int l4, m4, n4; // variable flag
  int l5, m5, n5; // variable function name
  int ierr = 0;

  doublereal out_f;
  logical    out_flag;

  GetRhsVar(1,"c",&m5, &n5, &l5); // Variable function name
  GetRhsVar(2,"i",&m1, &n1, &l1); // Variable n
  GetRhsVar(3,"d",&m2, &n2, &l2); // Variable x(n)

  if (strcmp(cstk(l5),"america")==0)
    {
      ierr = america_evalf__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     &out_f,     // Variable f
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"bratu3db")==0)
    {
      ierr = bratu3db_evalf__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"cache")==0)
    {
      ierr = cache_evalf__(istk(l1),   // Variable n
			   dstk(l2),   // Variable x(n)
			   &out_f,     // Variable f
			   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"condor")==0)
    {
      ierr = condor_evalf__(istk(l1),   // Variable n
			    dstk(l2),   // Variable x(n)
			    &out_f,     // Variable f
			    &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"contor2")==0)
    {
      ierr = contor2_evalf__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     &out_f,     // Variable f
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"contor")==0)
    {
      ierr = contor_evalf__(istk(l1),   // Variable n
			    dstk(l2),   // Variable x(n)
			    &out_f,     // Variable f
			    &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"ellipsoid")==0)
    {
      ierr = ellipsoid_evalf__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       &out_f,     // Variable f
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"genpack-cc-mina")==0)
    {
      ierr = genpack_cc_mina_evalf__(istk(l1),   // Variable n
			             dstk(l2),   // Variable x(n)
				     &out_f,     // Variable f
				     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"genpack-csq-mina")==0)
    {
      ierr = genpack_csq_mina_evalf__(istk(l1),   // Variable n
			              dstk(l2),   // Variable x(n)
				      &out_f,     // Variable f
				      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"hardcube")==0)
    {
      ierr = hardcube_evalf__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"hardspheres")==0)
    {
      ierr = hardspheres_evalf__(istk(l1),   // Variable n
			         dstk(l2),   // Variable x(n)
				 &out_f,     // Variable f
				 &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"kissing2")==0)
    {
      ierr = kissing2_evalf__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"kissing")==0)
    {
      ierr = kissing_evalf__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     &out_f,     // Variable f
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"location")==0)
    {
      ierr = location_evalf__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"mountain1")==0)
    {
      ierr = mountain1_evalf__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       &out_f,     // Variable f
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"mountain2")==0)
    {
      ierr = mountain2_evalf__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       &out_f,     // Variable f
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"packccmn")==0)
    {
      ierr = packccmn_evalf__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"packccmn-feas")==0)
    {
      ierr = packccmn_feas_evalf__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
				   &out_f,     // Variable f
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"packcrmn-feas")==0)
    {
      ierr = packcrmn_feas_evalf__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
				   &out_f,     // Variable f
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"pedazos4")==0)
    {
      ierr = pedazos4_evalf__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"piecefit")==0)
    {
      ierr = piecefit_evalf__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"simfock2")==0)
    {
      ierr = simfock2_evalf__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"simfock")==0)
    {
      ierr = simfock_evalf__(istk(l1),    // Variable n
			     dstk(l2),    // Variable x(n)
			      &out_f,     // Variable f
			      &out_flag); // Variable flag
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  m4 = 1; n4 = 1; // Variable f
  m5 = 1; n5 = 1; // Variable flag

  CreateVar(4, "d",&m4, &n4, &l4); // Variable f
  CreateVar(5, "i",&m5, &n5, &l5); // Variable flag

  *dstk(l4) = out_f;
  *istk(l5) = out_flag;

  LhsVar(1) = 3;
  LhsVar(2) = 4;

  PutLhsVar();

  return 0;
}

/////////////////////
// interface evalg //
/////////////////////

int interface_evalg_(char * fname)
{
  int l1, m1, n1; // variable n
  int l2, m2, n2; // variable x(n)
  int l3, m3, n3; // variable g(n)
  int l4, m4, n4; // variable flag
  int l5, m5, n5; // variable function name

  int ierr = 0, i;

  doublereal * out_g;
  integer      out_flag;

  GetRhsVar(1,"c",&m5, &n5, &l5); // Variable function name
  GetRhsVar(2,"i",&m1, &n1, &l1); // Variable n
  GetRhsVar(3,"d",&m2, &n2, &l2); // Variable x(n)

  out_g = (doublereal *)malloc(*istk(l1) * sizeof(doublereal));

  if (strcmp(cstk(l5),"america")==0)
    {
      ierr = america_evalg__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     out_g,      // Variable g(n)
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"bratu3db")==0)
    {
      ierr = bratu3db_evalg__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_g,      // Variable g(n)
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"cache")==0)
    {
      ierr = cache_evalg__(istk(l1),   // Variable n
			   dstk(l2),   // Variable x(n)
			   out_g,      // Variable g(n)
			   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"condor")==0)
    {
      ierr = condor_evalg__(istk(l1),   // Variable n
			    dstk(l2),   // Variable x(n)
			    out_g,      // Variable g(n)
			    &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"contor2")==0)
    {
      ierr = contor2_evalg__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     out_g,      // Variable g(n)
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"contor")==0)
    {
      ierr = contor_evalg__(istk(l1),   // Variable n
			    dstk(l2),   // Variable x(n)
			    out_g,      // Variable g(n)
			    &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"ellipsoid")==0)
    {
      ierr = ellipsoid_evalg__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       out_g,      // Variable g(n)
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"genpack-cc-mina")==0)
    {
      ierr = genpack_cc_mina_evalg__(istk(l1),   // Variable n
			             dstk(l2),   // Variable x(n)
				     out_g,      // Variable g(n)
				     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"genpack-csq-mina")==0)
    {
      ierr = genpack_csq_mina_evalg__(istk(l1),   // Variable n
			              dstk(l2),   // Variable x(n)
				      out_g,      // Variable g(n)
				      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"hardcube")==0)
    {
      ierr = hardcube_evalg__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_g,      // Variable g(n)
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"hardspheres")==0)
    {
      ierr = hardspheres_evalg__(istk(l1),   // Variable n
			         dstk(l2),   // Variable x(n)
				 out_g,      // Variable g(n)
				 &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"kissing2")==0)
    {
      ierr = kissing2_evalg__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_g,      // Variable g(n)
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"kissing")==0)
    {
      ierr = kissing_evalg__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     out_g,      // Variable g(n)
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"location")==0)
    {
      ierr = location_evalg__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_g,      // Variable g(n)
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"mountain1")==0)
    {
      ierr = mountain1_evalg__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       out_g,      // Variable g(n)
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"mountain2")==0)
    {
      ierr = mountain2_evalg__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       out_g,      // Variable g(n)
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"packccmn")==0)
    {
      ierr = packccmn_evalg__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_g,      // Variable g(n)
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"packccmn-feas")==0)
    {
      ierr = packccmn_feas_evalg__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
				   out_g,      // Variable g(n)
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"packcrmn-feas")==0)
    {
      ierr = packcrmn_feas_evalg__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
				   out_g,      // Variable g(n)
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"pedazos4")==0)
    {
      ierr = pedazos4_evalg__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_g,      // Variable g(n)
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"piecefit")==0)
    {
      ierr = piecefit_evalg__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_g,      // Variable g(n)
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"simfock2")==0)
    {
      ierr = simfock2_evalg__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_g,      // Variable g(n)
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l5),"simfock")==0)
    {
      ierr = simfock_evalg__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     out_g,      // Variable g(n)
			     &out_flag); // Variable flag
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  m4 = *istk(l1); n4 = 1; // Variable g(n)
  m5 = 1; n5 = 1; // Variable flag

  CreateVar(4, "d",&m4, &n4, &l4); // Variable g(n)
  CreateVar(5, "i",&m5, &n5, &l5); // Variable flag

  for(i=0; i<*istk(l1); i++)
    {
      *(dstk(l4)+i) = out_g[i];
    } /* End For */
  *istk(l5) = out_flag;

  LhsVar(1) = 3;
  LhsVar(2) = 4;

  PutLhsVar();

  return 0;
}

/////////////////////
// interface evalh //
/////////////////////

int interface_evalh_(char * fname)
{
  int l1, m1, n1; // variable n
  int l2, m2, n2; // variable x(n)
  int l3, m3, n3; // variable hnnz
  int l4, m4, n4; // variable hlin(hnnz)
  int l5, m5, n5; // variable hcol(hnnz)
  int l6, m6, n6; // variable hval(nnz)
  int l7, m7, n7; // variable flag
  int l8, m8, n8; // variable fname
  
  int ierr = 0, i;

  logical out_flag;
  integer out_nnzh;
  integer * out_hlin, * out_hcol;
  doublereal * out_hval;

  GetRhsVar(1,"c",&m8, &n8, &l8); // Variable function name
  GetRhsVar(2,"i",&m1, &n1, &l1); // Variable n
  GetRhsVar(3,"d",&m2, &n2, &l2); // Variable x(n)

  if (strcmp(cstk(l8),"america")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = america_evalh__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     out_hlin,   // Variable hlin
			     out_hcol,   // Variable hcol
			     out_hval,   // Variable hval
			     &out_nnzh,  // Variable nnzh
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"bratu3db")==0)
    {
      out_nnzh = 1000; // tmax in america.f
      out_hlin = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hcol = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hval = (doublereal *)malloc(out_nnzh * sizeof(doublereal));

      ierr = bratu3db_evalh__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_hlin,   // Variable hlin
			      out_hcol,   // Variable hcol
			      out_hval,   // Variable hval
			      &out_nnzh,  // Variable nnzh
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"cache")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = cache_evalh__(istk(l1),   // Variable n
			   dstk(l2),   // Variable x(n)
			   out_hlin,   // Variable hlin
			   out_hcol,   // Variable hcol
			   out_hval,   // Variable hval
			   &out_nnzh,  // Variable nnzh
			   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"condor")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = condor_evalh__(istk(l1),   // Variable n
			    dstk(l2),   // Variable x(n)
			    out_hlin,   // Variable hlin
			    out_hcol,   // Variable hcol
			    out_hval,   // Variable hval
			    &out_nnzh,  // Variable nnzh
			    &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"contor2")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = contor2_evalh__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     out_hlin,   // Variable hlin
			     out_hcol,   // Variable hcol
			     out_hval,   // Variable hval
			     &out_nnzh,  // Variable nnzh
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"contor")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = contor_evalh__(istk(l1),   // Variable n
			    dstk(l2),   // Variable x(n)
			    out_hlin,   // Variable hlin
			    out_hcol,   // Variable hcol
			    out_hval,   // Variable hval
			    &out_nnzh,  // Variable nnzh
			    &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"ellipsoid")==0)
    {
      out_nnzh = *istk(l1);
      out_hlin = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hcol = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hval = (doublereal *)malloc(out_nnzh * sizeof(doublereal));
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = ellipsoid_evalh__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       out_hlin,   // Variable hlin
			       out_hcol,   // Variable hcol
			       out_hval,   // Variable hval
			       &out_nnzh,  // Variable nnzh
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"genpack-cc-mina")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = genpack_cc_mina_evalh__(istk(l1),   // Variable n
			             dstk(l2),   // Variable x(n)
				     out_hlin,   // Variable hlin
				     out_hcol,   // Variable hcol
				     out_hval,   // Variable hval
				     &out_nnzh,  // Variable nnzh
				     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"genpack-csq-mina")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = genpack_csq_mina_evalh__(istk(l1),   // Variable n
			              dstk(l2),   // Variable x(n)
				      out_hlin,   // Variable hlin
				      out_hcol,   // Variable hcol
				      out_hval,   // Variable hval
				      &out_nnzh,  // Variable nnzh
				      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"hardcube")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = hardcube_evalh__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_hlin,   // Variable hlin
			      out_hcol,   // Variable hcol
			      out_hval,   // Variable hval
			      &out_nnzh,  // Variable nnzh
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"hardspheres")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = hardspheres_evalh__(istk(l1),   // Variable n
			         dstk(l2),   // Variable x(n)
				 out_hlin,   // Variable hlin
				 out_hcol,   // Variable hcol
				 out_hval,   // Variable hval
				 &out_nnzh,  // Variable nnzh
				 &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"kissing2")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = kissing2_evalh__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_hlin,   // Variable hlin
			      out_hcol,   // Variable hcol
			      out_hval,   // Variable hval
			      &out_nnzh,  // Variable nnzh
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"kissing")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = kissing_evalh__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     out_hlin,   // Variable hlin
			     out_hcol,   // Variable hcol
			     out_hval,   // Variable hval
			     &out_nnzh,  // Variable nnzh
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"location")==0)
    {
      out_nnzh = 10 * *istk(l1); // 10*n
      out_hlin = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hcol = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hval = (doublereal *)malloc(out_nnzh * sizeof(doublereal));

      ierr = location_evalh__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_hlin,   // Variable hlin
			      out_hcol,   // Variable hcol
			      out_hval,   // Variable hval
			      &out_nnzh,  // Variable nnzh
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"mountain1")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = mountain1_evalh__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       out_hlin,   // Variable hlin
			       out_hcol,   // Variable hcol
			       out_hval,   // Variable hval
			       &out_nnzh,  // Variable nnzh
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"mountain2")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = mountain2_evalh__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       out_hlin,   // Variable hlin
			       out_hcol,   // Variable hcol
			       out_hval,   // Variable hval
			       &out_nnzh,  // Variable nnzh
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"packccmn")==0)
    {
      out_nnzh = 80 * 2 * 2; // nite * ndim * ndim
      out_hlin = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hcol = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hval = (doublereal *)malloc(out_nnzh * sizeof(doublereal));

      ierr = packccmn_evalh__(istk(l1),    // Variable n
			      dstk(l2),    // Variable x(n)
			       out_hlin,   // Variable hlin
			       out_hcol,   // Variable hcol
			       out_hval,   // Variable hval
			       &out_nnzh,  // Variable nnzh
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"packccmn-feas")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = packccmn_feas_evalh__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
				   out_hlin,   // Variable hlin
				   out_hcol,   // Variable hcol
				   out_hval,   // Variable hval
				   &out_nnzh,  // Variable nnzh
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"packcrmn-feas")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = packcrmn_feas_evalh__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
				   out_hlin,   // Variable hlin
				   out_hcol,   // Variable hcol
				   out_hval,   // Variable hval
				   &out_nnzh,  // Variable nnzh
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"pedazos4")==0)
    {
      out_nnzh = 19; // nite * ndim * ndim
      out_hlin = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hcol = (integer *)malloc(out_nnzh * sizeof(integer));
      out_hval = (doublereal *)malloc(out_nnzh * sizeof(doublereal));

      ierr = pedazos4_evalh__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_hlin,   // Variable hlin
			      out_hcol,   // Variable hcol
			      out_hval,   // Variable hval
			      &out_nnzh,  // Variable nnzh
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"piecefit")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = piecefit_evalh__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_hlin,   // Variable hlin
			      out_hcol,   // Variable hcol
			      out_hval,   // Variable hval
			      &out_nnzh,  // Variable nnzh
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"simfock2")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = simfock2_evalh__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      out_hlin,   // Variable hlin
			      out_hcol,   // Variable hcol
			      out_hval,   // Variable hval
			      &out_nnzh,  // Variable nnzh
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l8),"simfock")==0)
    {
      out_hlin = (integer *)malloc(1 * sizeof(integer));
      out_hcol = (integer *)malloc(1 * sizeof(integer));
      out_hval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzh = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = simfock_evalh__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     out_hlin,   // Variable hlin
			     out_hcol,   // Variable hcol
			     out_hval,   // Variable hval
			     &out_nnzh,  // Variable nnzh
			     &out_flag); // Variable flag
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  m3 = 1; n3 = 1; // Variable hnnz
  m4 = out_nnzh; n4 = 1; // Variable hlin(hnnz)
  m5 = out_nnzh; n5 = 1; // Variable hcol(hnnz)
  m6 = out_nnzh; n6 = 1; // Variable hval(hnnz)
  m7 = 1; n7 = 1; // Variable flag

  CreateVar(4, "i",&m3, &n3, &l3); // Variable hnnz
  CreateVar(5, "i",&m4, &n4, &l4); // Variable hlin(hnnz)
  CreateVar(6, "i",&m5, &n5, &l5); // Variable hcol(hnnz)
  CreateVar(7, "d",&m6, &n6, &l6); // Variable hval(hnnz)
  CreateVar(8, "i",&m7, &n7, &l7); // Variable flag

  for(i=0; i<out_nnzh; i++)
    {
      *(istk(l4)+i) = out_hlin[i];
      *(istk(l5)+i) = out_hcol[i];
      *(dstk(l6)+i) = out_hval[i];
    } 
  *istk(l3) = out_nnzh;
  *istk(l7) = out_flag;

  free(out_hlin);
  free(out_hcol);
  free(out_hval);

  LhsVar(1) = 4;
  LhsVar(2) = 5;
  LhsVar(3) = 6;
  LhsVar(4) = 7;
  LhsVar(5) = 8;

  PutLhsVar();

  return 0;
}

/////////////////////
// interface evalc //
/////////////////////

int interface_evalc_(char * fname)
{
  int l1, m1, n1; // variable n
  int l2, m2, n2; // variable x(n)
  int l3, m3, n3; // variable ind
  int l4, m4, n4; // variable c
  int l5, m5, n5; // variable flag
  int l6, m6, n6; // variable fname
  
  int ierr = 0;

  doublereal out_c;
  logical    out_flag;
  integer    in_index;

  GetRhsVar(1,"c",&m6, &n6, &l6); // Variable function name
  GetRhsVar(2,"i",&m1, &n1, &l1); // Variable n
  GetRhsVar(3,"d",&m2, &n2, &l2); // Variable x(n)
  GetRhsVar(4,"i",&m3, &n3, &l3); // Variable ind

  in_index = *(istk(l3));

  if (strcmp(cstk(l6),"america")==0)
    {
      ierr = america_evalc__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     &in_index,  // Variable ind
			     &out_c,     // Variable c
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"bratu3db")==0)
    {
      ierr = bratu3db_evalc__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable ind
			      &out_c,     // Variable c
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"cache")==0)
    {
      ierr = cache_evalc__(istk(l1),   // Variable n
			   dstk(l2),   // Variable x(n)
			   istk(l3),   // Variable ind
			   &out_c,     // Variable c
			   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"condor")==0)
    {
      ierr = condor_evalc__(istk(l1),   // Variable n
			    dstk(l2),   // Variable x(n)
			    istk(l3),   // Variable ind
			    &out_c,     // Variable c
			    &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"contor2")==0)
    {
      ierr = contor2_evalc__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     istk(l3),   // Variable ind
			     &out_c,     // Variable c
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"contor")==0)
    {
      ierr = condor_evalc__(istk(l1),   // Variable n
			    dstk(l2),   // Variable x(n)
			    istk(l3),   // Variable ind
			    &out_c,     // Variable c
			    &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"ellipsoid")==0)
    {
      ierr = ellipsoid_evalc__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       istk(l3),   // Variable ind
			       &out_c,     // Variable c
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"genpack-cc-mina")==0)
    {
      ierr = genpack_cc_mina_evalc__(istk(l1),   // Variable n
			             dstk(l2),   // Variable x(n)
			             istk(l3),   // Variable ind
				     &out_c,     // Variable c
				     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"genpack-csq-mina")==0)
    {
      ierr = genpack_csq_mina_evalc__(istk(l1),   // Variable n
		    	              dstk(l2),   // Variable x(n)
			              istk(l3),   // Variable ind
				      &out_c,     // Variable c
				      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"hardcube")==0)
    {
      ierr = hardcube_evalc__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable ind
			      &out_c,     // Variable c
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"hardspheres")==0)
    {
      ierr = hardspheres_evalc__(istk(l1),   // Variable n
			         dstk(l2),   // Variable x(n)
			         istk(l3),   // Variable ind
				 &out_c,     // Variable c
				 &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"kissing2")==0)
    {
      ierr = kissing2_evalc__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable ind
			      &out_c,     // Variable c
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"kissing")==0)
    {
      ierr = kissing_evalc__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     istk(l3),   // Variable ind
			     &out_c,     // Variable c
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"location")==0)
    {
      ierr = location_evalc__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable ind
			      &out_c,     // Variable c
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"mountain1")==0)
    {
      ierr = mountain1_evalc__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       istk(l3),   // Variable ind
			       &out_c,     // Variable c
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"mountain2")==0)
    {
      ierr = mountain2_evalc__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       istk(l3),   // Variable ind
			       &out_c,     // Variable c
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"packccmn")==0)
    {
      ierr = packccmn_evalc__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable ind
			      &out_c,     // Variable c
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"packccmn-feas")==0)
    {
      ierr = packccmn_feas_evalc__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
			           istk(l3),   // Variable ind
				   &out_c,     // Variable c
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"packcrmn-feas")==0)
    {
      ierr = packcrmn_feas_evalc__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
			           istk(l3),   // Variable ind
				   &out_c,     // Variable c
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"pedazos4")==0)
    {
      ierr = pedazos4_evalc__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable ind
			      &out_c,     // Variable c
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"piecefit")==0)
    {
      ierr = piecefit_evalc__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable ind
			      &out_c,     // Variable c
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"simfock2")==0)
    {
      ierr = simfock2_evalc__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable ind
			      &out_c,     // Variable c
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l6),"simfock")==0)
    {
      ierr = simfock_evalc__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     istk(l3),   // Variable ind
			     &out_c,     // Variable c
			     &out_flag); // Variable flag
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  m4 = 1; n4 = 1; // Variable c
  m5 = 1; n5 = 1; // Variable flag

  CreateVar(5, "d",&m4, &n4, &l4); // Variable c
  CreateVar(6, "i",&m5, &n5, &l5); // Variable flag

  *(dstk(l4)) = out_c;
  *(lstk(l5))    = out_flag;

  LhsVar(1) = 5;
  LhsVar(2) = 6;

  PutLhsVar();

  return 0;
}

///////////////////////
// interface evaljac //
///////////////////////

int interface_evaljac_(char * fname)
{
  int l1, m1, n1; // variable n
  int l2, m2, n2; // variable x(n)
  int l3, m3, n3; // variable ind
  int l4, m4, n4; // variable jcnnz
  int l5, m5, n5; // variable jcvar(jcnnz)
  int l6, m6, n6; // variable jcval(jcnnz)
  int l7, m7, n7; // variable flag
  int l8, m8, n8; // variable fname

  int ierr = 0, i;

  integer out_nnzjac;
  integer * out_indjac;
  doublereal * out_valjac;
  logical out_flag;

  GetRhsVar(1,"c",&m8, &n8, &l8); // Variable function name
  GetRhsVar(2,"i",&m1, &n1, &l1); // Variable n
  GetRhsVar(3,"d",&m2, &n2, &l2); // Variable x(n)
  GetRhsVar(4,"i",&m3, &n3, &l3); // Variable ind

  if (strcmp(cstk(l8),"america")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = america_evaljac__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_indjac,  // Variable indjac(nnzjac)
			       out_valjac,  // Variable valjac(nnzjac)
			       &out_nnzjac, // Variable nnzjac
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"bratu3db")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = bratu3db_evaljac__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_indjac,  // Variable indjac(nnzjac)
				out_valjac,  // Variable valjac(nnzjac)
				&out_nnzjac, // Variable nnzjac
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"cache")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = cache_evaljac__(istk(l1),    // Variable n
			     dstk(l2),    // Variable x(n)
			     istk(l3),    // Variable ind
			     out_indjac,  // Variable indjac(nnzjac)
			     out_valjac,  // Variable valjac(nnzjac)
			     &out_nnzjac, // Variable nnzjac
			     &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"condor")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = condor_evaljac__(istk(l1),    // Variable n
			      dstk(l2),    // Variable x(n)
			      istk(l3),    // Variable ind
			      out_indjac,  // Variable indjac(nnzjac)
			      out_valjac,  // Variable valjac(nnzjac)
			      &out_nnzjac, // Variable nnzjac
			      &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"contor2")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = contor2_evaljac__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_indjac,  // Variable indjac(nnzjac)
			       out_valjac,  // Variable valjac(nnzjac)
			       &out_nnzjac, // Variable nnzjac
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"contor")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = contor_evaljac__(istk(l1),    // Variable n
			      dstk(l2),    // Variable x(n)
			      istk(l3),    // Variable ind
			      out_indjac,  // Variable indjac(nnzjac)
			      out_valjac,  // Variable valjac(nnzjac)
			      &out_nnzjac, // Variable nnzjac
			      &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"ellipsoid")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = ellipsoid_evaljac__(istk(l1),    // Variable n
			         dstk(l2),    // Variable x(n)
			         istk(l3),    // Variable ind
				 out_indjac,  // Variable indjac(nnzjac)
				 out_valjac,  // Variable valjac(nnzjac)
				 &out_nnzjac, // Variable nnzjac
				 &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"genpack-cc-mina")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = genpack_cc_mina_evaljac__(istk(l1),    // Variable n
			               dstk(l2),    // Variable x(n)
			               istk(l3),    // Variable ind
				       out_indjac,  // Variable indjac(nnzjac)
				       out_valjac,  // Variable valjac(nnzjac)
				       &out_nnzjac, // Variable nnzjac
				       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"genpack-csq-mina")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = genpack_csq_mina_evaljac__(istk(l1),    // Variable n
			                dstk(l2),    // Variable x(n)
			                istk(l3),    // Variable ind
					out_indjac,  // Variable indjac(nnzjac)
					out_valjac,  // Variable valjac(nnzjac)
					&out_nnzjac, // Variable nnzjac
					&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"hardcube")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = hardcube_evaljac__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_indjac,  // Variable indjac(nnzjac)
				out_valjac,  // Variable valjac(nnzjac)
				&out_nnzjac, // Variable nnzjac
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"hardspheres")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = hardspheres_evaljac__(istk(l1),    // Variable n
			           dstk(l2),    // Variable x(n)
			           istk(l3),    // Variable ind
				   out_indjac,  // Variable indjac(nnzjac)
				   out_valjac,  // Variable valjac(nnzjac)
				   &out_nnzjac, // Variable nnzjac
				   &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"kissing2")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = kissing2_evaljac__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_indjac,  // Variable indjac(nnzjac)
				out_valjac,  // Variable valjac(nnzjac)
				&out_nnzjac, // Variable nnzjac
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"kissing")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = kissing_evaljac__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_indjac,  // Variable indjac(nnzjac)
			       out_valjac,  // Variable valjac(nnzjac)
			       &out_nnzjac, // Variable nnzjac
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"location")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = location_evaljac__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_indjac,  // Variable indjac(nnzjac)
				out_valjac,  // Variable valjac(nnzjac)
				&out_nnzjac, // Variable nnzjac
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"mountain1")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = mountain1_evaljac__(istk(l1),    // Variable n
			         dstk(l2),    // Variable x(n)
			         istk(l3),    // Variable ind
				 out_indjac,  // Variable indjac(nnzjac)
				 out_valjac,  // Variable valjac(nnzjac)
				 &out_nnzjac, // Variable nnzjac
				 &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"mountain2")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = mountain2_evaljac__(istk(l1),    // Variable n
			         dstk(l2),    // Variable x(n)
			         istk(l3),    // Variable ind
				 out_indjac,  // Variable indjac(nnzjac)
				 out_valjac,  // Variable valjac(nnzjac)
				 &out_nnzjac, // Variable nnzjac
				 &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"packccmn")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = packccmn_evaljac__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_indjac,  // Variable indjac(nnzjac)
				out_valjac,  // Variable valjac(nnzjac)
				&out_nnzjac, // Variable nnzjac
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"packccmn-feas")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = packccmn_feas_evaljac__(istk(l1),    // Variable n
			             dstk(l2),    // Variable x(n)
			             istk(l3),    // Variable ind
				     out_indjac,  // Variable indjac(nnzjac)
				     out_valjac,  // Variable valjac(nnzjac)
				     &out_nnzjac, // Variable nnzjac
				     &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"packcrmn-feas")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = packcrmn_feas_evaljac__(istk(l1),    // Variable n
			             dstk(l2),    // Variable x(n)
			             istk(l3),    // Variable ind
				     out_indjac,  // Variable indjac(nnzjac)
				     out_valjac,  // Variable valjac(nnzjac)
				     &out_nnzjac, // Variable nnzjac
				     &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"pedazos4")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = pedazos4_evaljac__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_indjac,  // Variable indjac(nnzjac)
				out_valjac,  // Variable valjac(nnzjac)
				&out_nnzjac, // Variable nnzjac
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"piecefit")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = piecefit_evaljac__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_indjac,  // Variable indjac(nnzjac)
				out_valjac,  // Variable valjac(nnzjac)
				&out_nnzjac, // Variable nnzjac
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"simfock2")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = simfock2_evaljac__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_indjac,  // Variable indjac(nnzjac)
				out_valjac,  // Variable valjac(nnzjac)
				&out_nnzjac, // Variable nnzjac
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"simfock")==0)
    {
      out_nnzjac = *istk(l1);
      out_indjac = (integer *)malloc(out_nnzjac * sizeof(integer));
      out_valjac = (doublereal *)malloc(out_nnzjac * sizeof(doublereal));

      ierr = simfock_evaljac__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_indjac,  // Variable indjac(nnzjac)
			       out_valjac,  // Variable valjac(nnzjac)
			       &out_nnzjac, // Variable nnzjac
			       &out_flag);  // Variable flag
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  m4 = 1; n4 = 1;   // Variable nnzjac
  m5 = out_nnzjac; n5 = 1; // Variable indjac(nnzjac)
  m6 = out_nnzjac; n6 = 1; // Variable valjac(nnzjac)
  m7 = 1; n7 = 1;   // Variable flag

  CreateVar(5, "i",&m4, &n4, &l4); // Variable nnzjac
  CreateVar(6, "i",&m5, &n5, &l5); // Variable indjac(jcnnz)
  CreateVar(7, "d",&m6, &n6, &l6); // Variable valjac(jcnnz)
  CreateVar(8, "i",&m7, &n7, &l7); // Variable flag(jcnnz)

  for(i=0; i<out_nnzjac; i++)
    {
      *(istk(l5)+i) = out_indjac[i];
      *(dstk(l6)+i) = out_valjac[i];
    } 
  *istk(l4) = out_nnzjac;
  *istk(l7) = out_flag;

  free(out_indjac);
  free(out_valjac);

  LhsVar(1) = 5;
  LhsVar(2) = 6;
  LhsVar(3) = 7;
  LhsVar(4) = 8;

  PutLhsVar();

  return 0;
}

//////////////////////
// interface evalhc //
//////////////////////

int interface_evalhc_(char * fname)
{
  int l1, m1, n1; // variable n
  int l2, m2, n2; // variable x(n)
  int l3, m3, n3; // variable ind
  int l4, m4, n4; // variable hcnnz
  int l5, m5, n5; // variable hclin(jcnnz)
  int l6, m6, n6; // variable hccol(jcnnz)
  int l7, m7, n7; // variable hcval(hcnnz)
  int l8, m8, n8; // variable flag
  int l9, m9, n9; // variable fname

  int ierr = 0, i;

  integer out_nnzhc;
  logical out_flag;
  integer * out_hclin, * out_hccol;
  doublereal * out_hcval;

  GetRhsVar(1,"c",&m9, &n9, &l9); // Variable function name
  GetRhsVar(2,"i",&m1, &n1, &l1); // Variable n
  GetRhsVar(3,"d",&m2, &n2, &l2); // Variable x(n)
  GetRhsVar(4,"i",&m3, &n3, &l3); // Variable ind

  if (strcmp(cstk(l9),"america")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = america_evalhc__(istk(l1),    // Variable n
			      dstk(l2),    // Variable x(n)
			      istk(l3),    // Variable ind
			      out_hclin,   // Variable hclin(nnzhc)
			      out_hccol,   // Variable hccol(nnzhc)
			      out_hcval,   // Variable hcval(nnzhc)
			      &out_nnzhc,  // Variable nnzhc
			      &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"bratu3db")==0)
    {
      out_nnzhc = 1;
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = bratu3db_evalhc__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_hclin,   // Variable hclin(nnzhc)
			       out_hccol,   // Variable hccol(nnzhc)
			       out_hcval,   // Variable hcval(nnzhc)
			       &out_nnzhc,  // Variable nnzhc
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"cache")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = cache_evalhc__(istk(l1),    // Variable n
			    dstk(l2),    // Variable x(n)
			    istk(l3),    // Variable ind
			    out_hclin,   // Variable hclin(nnzhc)
			    out_hccol,   // Variable hccol(nnzhc)
			    out_hcval,   // Variable hcval(nnzhc)
			    &out_nnzhc,  // Variable nnzhc
			    &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"condor")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = condor_evalhc__(istk(l1),    // Variable n
			     dstk(l2),    // Variable x(n)
			     istk(l3),    // Variable ind
			     out_hclin,   // Variable hclin(nnzhc)
			     out_hccol,   // Variable hccol(nnzhc)
			     out_hcval,   // Variable hcval(nnzhc)
			     &out_nnzhc,  // Variable nnzhc
			     &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"contor2")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = contor2_evalhc__(istk(l1),    // Variable n
			      dstk(l2),    // Variable x(n)
			      istk(l3),    // Variable ind
			      out_hclin,   // Variable hclin(nnzhc)
			      out_hccol,   // Variable hccol(nnzhc)
			      out_hcval,   // Variable hcval(nnzhc)
			      &out_nnzhc,  // Variable nnzhc
			      &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"contor")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = contor_evalhc__(istk(l1),    // Variable n
			     dstk(l2),    // Variable x(n)
			     istk(l3),    // Variable ind
			     out_hclin,   // Variable hclin(nnzhc)
			     out_hccol,   // Variable hccol(nnzhc)
			     out_hcval,   // Variable hcval(nnzhc)
			     &out_nnzhc,  // Variable nnzhc
			     &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"ellipsoid")==0)
    {
      out_nnzhc = *istk(l1) * *istk(l1) * *istk(l1); // n^3
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = ellipsoid_evalhc__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_hclin,   // Variable hclin(nnzhc)
				out_hccol,   // Variable hccol(nnzhc)
				out_hcval,   // Variable hcval(nnzhc)
				&out_nnzhc,  // Variable nnzhc
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"genpack-cc-mina")==0)
    {
      out_nnzhc = *istk(l1) * *istk(l1) * *istk(l1); // n^3
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = genpack_cc_mina_evalhc__(istk(l1),    // Variable n
			              dstk(l2),    // Variable x(n)
			              istk(l3),    // Variable ind
				      out_hclin,   // Variable hclin(nnzhc)
				      out_hccol,   // Variable hccol(nnzhc)
				      out_hcval,   // Variable hcval(nnzhc)
				      &out_nnzhc,  // Variable nnzhc
				      &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"genpack-csq-mina")==0)
    {
      out_nnzhc = *istk(l1) * *istk(l1) * *istk(l1); // n^3
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = genpack_csq_mina_evalhc__(istk(l1),    // Variable n
			               dstk(l2),    // Variable x(n)
			               istk(l3),    // Variable ind
				       out_hclin,   // Variable hclin(nnzhc)
				       out_hccol,   // Variable hccol(nnzhc)
				       out_hcval,   // Variable hcval(nnzhc)
				       &out_nnzhc,  // Variable nnzhc
				       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"hardcube")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = hardcube_evalhc__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_hclin,   // Variable hclin(nnzhc)
			       out_hccol,   // Variable hccol(nnzhc)
			       out_hcval,   // Variable hcval(nnzhc)
			       &out_nnzhc,  // Variable nnzhc
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"hardspheres")==0)
    {
      out_nnzhc = *istk(l1); // n
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = hardspheres_evalhc__(istk(l1),    // Variable n
			          dstk(l2),    // Variable x(n)
			          istk(l3),    // Variable ind
				  out_hclin,   // Variable hclin(nnzhc)
				  out_hccol,   // Variable hccol(nnzhc)
				  out_hcval,   // Variable hcval(nnzhc)
				  &out_nnzhc,  // Variable nnzhc
				  &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"kissing2")==0)
    {
      out_nnzhc = *istk(l1); // n
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = kissing2_evalhc__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_hclin,   // Variable hclin(nnzhc)
			       out_hccol,   // Variable hccol(nnzhc)
			       out_hcval,   // Variable hcval(nnzhc)
			       &out_nnzhc,  // Variable nnzhc
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"kissing")==0)
    {
      out_nnzhc = *istk(l1); // n
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = kissing_evalhc__(istk(l1),    // Variable n
			      dstk(l2),    // Variable x(n)
			      istk(l3),    // Variable ind
			      out_hclin,   // Variable hclin(nnzhc)
			      out_hccol,   // Variable hccol(nnzhc)
			      out_hcval,   // Variable hcval(nnzhc)
			      &out_nnzhc,  // Variable nnzhc
			      &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"location")==0)
    {
      out_nnzhc = 2; 
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = location_evalhc__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_hclin,   // Variable hclin(nnzhc)
			       out_hccol,   // Variable hccol(nnzhc)
			       out_hcval,   // Variable hcval(nnzhc)
			       &out_nnzhc,  // Variable nnzhc
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"mountain1")==0)
    {
      out_nnzhc = 6; 
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = mountain1_evalhc__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_hclin,   // Variable hclin(nnzhc)
				out_hccol,   // Variable hccol(nnzhc)
				out_hcval,   // Variable hcval(nnzhc)
				&out_nnzhc,  // Variable nnzhc
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"mountain2")==0)
    {
      out_nnzhc = 6; 
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = mountain2_evalhc__(istk(l1),    // Variable n
			        dstk(l2),    // Variable x(n)
			        istk(l3),    // Variable ind
				out_hclin,   // Variable hclin(nnzhc)
				out_hccol,   // Variable hccol(nnzhc)
				out_hcval,   // Variable hcval(nnzhc)
				&out_nnzhc,  // Variable nnzhc
				&out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"packccmn")==0)
    {
      out_nnzhc = *istk(l1); 
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = packccmn_evalhc__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_hclin,   // Variable hclin(nnzhc)
			       out_hccol,   // Variable hccol(nnzhc)
			       out_hcval,   // Variable hcval(nnzhc)
			       &out_nnzhc,  // Variable nnzhc
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"packccmn-feas")==0)
    {
      out_nnzhc = *istk(l1); 
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = packccmn_feas_evalhc__(istk(l1),    // Variable n
			            dstk(l2),    // Variable x(n)
			            istk(l3),    // Variable ind
				    out_hclin,   // Variable hclin(nnzhc)
				    out_hccol,   // Variable hccol(nnzhc)
				    out_hcval,   // Variable hcval(nnzhc)
				    &out_nnzhc,  // Variable nnzhc
				    &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"packcrmn-feas")==0)
    {
      out_nnzhc = *istk(l1) * 3; 
      out_hclin = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hccol = (integer *)malloc(out_nnzhc * sizeof(integer));
      out_hcval = (doublereal *)malloc(out_nnzhc * sizeof(doublereal));

      ierr = packcrmn_feas_evalhc__(istk(l1),    // Variable n
			            dstk(l2),    // Variable x(n)
	     		            istk(l3),    // Variable ind
				    out_hclin,   // Variable hclin(nnzhc)
				    out_hccol,   // Variable hccol(nnzhc)
				    out_hcval,   // Variable hcval(nnzhc)
				    &out_nnzhc,  // Variable nnzhc
				    &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"pedazos4")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = pedazos4_evalhc__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_hclin,   // Variable hclin(nnzhc)
			       out_hccol,   // Variable hccol(nnzhc)
			       out_hcval,   // Variable hcval(nnzhc)
			       &out_nnzhc,  // Variable nnzhc
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"piecefit")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = piecefit_evalhc__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_hclin,   // Variable hclin(nnzhc)
			       out_hccol,   // Variable hccol(nnzhc)
			       out_hcval,   // Variable hcval(nnzhc)
			       &out_nnzhc,  // Variable nnzhc
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l9),"simfock2")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = simfock2_evalhc__(istk(l1),    // Variable n
			       dstk(l2),    // Variable x(n)
			       istk(l3),    // Variable ind
			       out_hclin,   // Variable hclin(nnzhc)
			       out_hccol,   // Variable hccol(nnzhc)
			       out_hcval,   // Variable hcval(nnzhc)
			       &out_nnzhc,  // Variable nnzhc
			       &out_flag);  // Variable flag
    }
  else if (strcmp(cstk(l8),"simfock")==0)
    {
      out_hclin = (integer *)malloc(1 * sizeof(integer));
      out_hccol = (integer *)malloc(1 * sizeof(integer));
      out_hcval = (doublereal *)malloc(1 * sizeof(doublereal));
      out_nnzhc = 0;
      // flag will be equal to -1: this function is not implemented for 'america'

      ierr = simfock_evalhc__(istk(l1),    // Variable n
			      dstk(l2),    // Variable x(n)
			      istk(l3),    // Variable ind
			      out_hclin,   // Variable hclin(nnzhc)
			      out_hccol,   // Variable hccol(nnzhc)
			      out_hcval,   // Variable hcval(nnzhc)
			      &out_nnzhc,  // Variable nnzhc
			      &out_flag);  // Variable flag
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  m4 = 1;  n4 = 1; // Variable hcnnz
  m5 = out_nnzhc; n5 = 1; // Variable hclin(hcnnz)
  m6 = out_nnzhc; n6 = 1; // Variable hccol(hcnnz)
  m7 = out_nnzhc; n7 = 1; // Variable hcvol(hcnnz)
  m8 = 1;  n8 = 1; // Variable flag

  CreateVar(5, "i",&m4, &n4, &l4); // Variable hcnnz
  CreateVar(6, "i",&m5, &n5, &l5); // Variable hclin(jcnnz)
  CreateVar(7, "i",&m6, &n6, &l6); // Variable hccol(jcnnz)
  CreateVar(8, "d",&m7, &n7, &l7); // Variable hcval(jcnnz)
  CreateVar(9, "i",&m8, &n8, &l8); // Variable flag

  for(i=0; i<out_nnzhc; i++)
    {
      *(istk(l5)+i) = out_hclin[i];
      *(istk(l6)+i) = out_hccol[i];
      *(dstk(l7)+i) = out_hcval[i];
    } 
  *istk(l4) = out_nnzhc;
  *istk(l8) = out_flag;

  free(out_hclin);
  free(out_hccol);
  free(out_hcval);

  LhsVar(1) = 5;
  LhsVar(2) = 6;
  LhsVar(3) = 7;
  LhsVar(4) = 8;
  LhsVar(5) = 9;

  PutLhsVar();

  return 0;
}

///////////////////////
// interface evalhlp //
///////////////////////

int interface_evalhlp_(char * fname)
{
  int l1, m1, n1; // variable n
  int l2, m2, n2; // variable x(n)
  int l3, m3, n3; // variable m
  int l4, m4, n4; // variable lambda(m)
  int l5, m5, n5; // variable p(n)
  int l6, m6, n6; // variable goth
  int l7, m7, n7; // variable hp(n)
  int l8, m8, n8; // variable flag
  int l9, m9, n9; // variable fname

  int ierr = 0, i;

  doublereal * out_hp;
  logical      out_flag;

  GetRhsVar(1,"c",&m9, &n9, &l9); // Variable function name
  GetRhsVar(2,"i",&m1, &n1, &l1); // Variable n
  GetRhsVar(3,"d",&m2, &n2, &l2); // Variable x(n)
  GetRhsVar(4,"i",&m3, &n3, &l3); // Variable m
  GetRhsVar(5,"d",&m4, &n4, &l4); // Variable lambda(m)
  GetRhsVar(6,"d",&m5, &n5, &l5); // Variable p(n)
  GetRhsVar(7,"i",&m6, &n6, &l6); // Variable goth

  out_hp = (doublereal *)malloc(*istk(l1) * sizeof(doublereal));

  if (strcmp(cstk(l9),"america")==0)
    {
      ierr = america_evalhlp__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       istk(l3),   // Variable m
			       dstk(l4),   // Variable lambda(m)
			       dstk(l5),   // Variable p(m)
			       out_hp,     // Variable hp(n)
			       lstk(l6),   // Variable goth
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"bratu3db")==0)
    {
      ierr = bratu3db_evalhlp__(istk(l1),   // Variable n
			        dstk(l2),   // Variable x(n)
			        istk(l3),   // Variable m
			        dstk(l4),   // Variable lambda(m)
			        dstk(l5),   // Variable p(m)
				out_hp,     // Variable hp(n)
				lstk(l6),   // Variable goth
				&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"cache")==0)
    {
      ierr = cache_evalhlp__(istk(l1),   // Variable n
			     dstk(l2),   // Variable x(n)
			     istk(l3),   // Variable m
			     dstk(l4),   // Variable lambda(m)
			     dstk(l5),   // Variable p(m)
			     out_hp,     // Variable hp(n)
			     lstk(l6),   // Variable goth
			     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"condor")==0)
    {
      ierr = condor_evalhlp__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable m
			      dstk(l4),   // Variable lambda(m)
			      dstk(l5),   // Variable p(m)
			      out_hp,     // Variable hp(n)
			      lstk(l6),   // Variable goth
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"contor2")==0)
    {
      ierr = contor2_evalhlp__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       istk(l3),   // Variable m
			       dstk(l4),   // Variable lambda(m)
			       dstk(l5),   // Variable p(m)
			       out_hp,     // Variable hp(n)
			       lstk(l6),   // Variable goth
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"contor")==0)
    {
      ierr = contor_evalhlp__(istk(l1),   // Variable n
			      dstk(l2),   // Variable x(n)
			      istk(l3),   // Variable m
			      dstk(l4),   // Variable lambda(m)
			      dstk(l5),   // Variable p(m)
			      out_hp,     // Variable hp(n)
			      lstk(l6),   // Variable goth
			      &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"ellipsoid")==0)
    {
      ierr = ellipsoid_evalhlp__(istk(l1),   // Variable n
			         dstk(l2),   // Variable x(n)
			         istk(l3),   // Variable m
			         dstk(l4),   // Variable lambda(m)
			         dstk(l5),   // Variable p(m)
				 out_hp,     // Variable hp(n)
				 lstk(l6),   // Variable goth
				 &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"genpack-cc-mina")==0)
    {
      ierr = genpack_cc_mina_evalhlp__(istk(l1),   // Variable n
			               dstk(l2),   // Variable x(n)
			               istk(l3),   // Variable m
			               dstk(l4),   // Variable lambda(m)
			               dstk(l5),   // Variable p(m)
				       out_hp,     // Variable hp(n)
				       lstk(l6),   // Variable goth
				       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"genpack-csq-mina")==0)
    {
      ierr = genpack_csq_mina_evalhlp__(istk(l1),   // Variable n
			                dstk(l2),   // Variable x(n)
			                istk(l3),   // Variable m
			                dstk(l4),   // Variable lambda(m)
			                dstk(l5),   // Variable p(m)
					out_hp,     // Variable hp(n)
					lstk(l6),   // Variable goth
					&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"hardcube")==0)
    {
      ierr = hardcube_evalhlp__(istk(l1),   // Variable n
			        dstk(l2),   // Variable x(n)
			        istk(l3),   // Variable m
			        dstk(l4),   // Variable lambda(m)
			        dstk(l5),   // Variable p(m)
				out_hp,     // Variable hp(n)
				lstk(l6),   // Variable goth
				&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"hardspheres")==0)
    {
      ierr = hardspheres_evalhlp__(istk(l1),   // Variable n
			           dstk(l2),   // Variable x(n)
			           istk(l3),   // Variable m
			           dstk(l4),   // Variable lambda(m)
			           dstk(l5),   // Variable p(m)
				   out_hp,     // Variable hp(n)
				   lstk(l6),   // Variable goth
				   &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"kissing2")==0)
    {
      ierr = kissing2_evalhlp__(istk(l1),   // Variable n
			        dstk(l2),   // Variable x(n)
			        istk(l3),   // Variable m
			        dstk(l4),   // Variable lambda(m)
			        dstk(l5),   // Variable p(m)
				out_hp,     // Variable hp(n)
				lstk(l6),   // Variable goth
				&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"kissing")==0)
    {
      ierr = kissing_evalhlp__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       istk(l3),   // Variable m
			       dstk(l4),   // Variable lambda(m)
			       dstk(l5),   // Variable p(m)
			       out_hp,     // Variable hp(n)
			       lstk(l6),   // Variable goth
			       &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"location")==0)
    {
      ierr = location_evalhlp__(istk(l1),   // Variable n
			        dstk(l2),   // Variable x(n)
			        istk(l3),   // Variable m
			        dstk(l4),   // Variable lambda(m)
			        dstk(l5),   // Variable p(m)
				out_hp,     // Variable hp(n)
				lstk(l6),   // Variable goth
				&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"mountain1")==0)
    {
      ierr = mountain1_evalhlp__(istk(l1),   // Variable n
				 dstk(l2),   // Variable x(n)
				 istk(l3),   // Variable m
				 dstk(l4),   // Variable lambda(m)
				 dstk(l5),   // Variable p(m)
				 out_hp,     // Variable hp(n)
				 lstk(l6),   // Variable goth
				 &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"mountain2")==0)
    {
      ierr = mountain2_evalhlp__(istk(l1),   // Variable n
				 dstk(l2),   // Variable x(n)
				 istk(l3),   // Variable m
				 dstk(l4),   // Variable lambda(m)
				 dstk(l5),   // Variable p(m)
				 out_hp,     // Variable hp(n)
				 lstk(l6),   // Variable goth
				 &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"packccmn")==0)
    {
      ierr = packccmn_evalhlp__(istk(l1),   // Variable n
			        dstk(l2),   // Variable x(n)
			        istk(l3),   // Variable m
			        dstk(l4),   // Variable lambda(m)
			        dstk(l5),   // Variable p(m)
				out_hp,     // Variable hp(n)
				lstk(l6),   // Variable goth
				&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"packccmn-feas")==0)
    {
      ierr = packccmn_feas_evalhlp__(istk(l1),   // Variable n
				     dstk(l2),   // Variable x(n)
				     istk(l3),   // Variable m
				     dstk(l4),   // Variable lambda(m)
				     dstk(l5),   // Variable p(m)
				     out_hp,     // Variable hp(n)
				     lstk(l6),   // Variable goth
				     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"packcrmn-feas")==0)
    {
      ierr = packcrmn_feas_evalhlp__(istk(l1),   // Variable n
				     dstk(l2),   // Variable x(n)
				     istk(l3),   // Variable m
				     dstk(l4),   // Variable lambda(m)
				     dstk(l5),   // Variable p(m)
				     out_hp,     // Variable hp(n)
				     lstk(l6),   // Variable goth
				     &out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"pedazos4")==0)
    {
      ierr = pedazos4_evalhlp__(istk(l1),   // Variable n
			        dstk(l2),   // Variable x(n)
			        istk(l3),   // Variable m
			        dstk(l4),   // Variable lambda(m)
			        dstk(l5),   // Variable p(m)
				out_hp,     // Variable hp(n)
				lstk(l6),   // Variable goth
				&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"piecefit")==0)
    {
      ierr = piecefit_evalhlp__(istk(l1),   // Variable n
			        dstk(l2),   // Variable x(n)
			        istk(l3),   // Variable m
			        dstk(l4),   // Variable lambda(m)
			        dstk(l5),   // Variable p(m)
				out_hp,     // Variable hp(n)
				lstk(l6),   // Variable goth
				&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"simfock2")==0)
    {
      ierr = simfock2_evalhlp__(istk(l1),   // Variable n
			        dstk(l2),   // Variable x(n)
			        istk(l3),   // Variable m
			        dstk(l4),   // Variable lambda(m)
			        dstk(l5),   // Variable p(m)
				out_hp,     // Variable hp(n)
				lstk(l6),   // Variable goth
				&out_flag); // Variable flag
    }
  else if (strcmp(cstk(l9),"simfock")==0)
    {
      ierr = simfock_evalhlp__(istk(l1),   // Variable n
			       dstk(l2),   // Variable x(n)
			       istk(l3),   // Variable m
			       dstk(l4),   // Variable lambda(m)
			       dstk(l5),   // Variable p(m)
			       out_hp,     // Variable hp(n)
			       lstk(l6),   // Variable goth
			       &out_flag); // Variable flag
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  m7 = *istk(l1); n7 = 1; // Variable hp(n)
  m6 = 1; n6 = 1; // Variable goth
  m8 = 1; n8 = 1; // Variable flag

  CreateVar(8, "d",&m7, &n7, &l7); // Variable hp(n)
  CreateVar(9, "i",&m6, &n6, &l6); // Variable goth
  CreateVar(10,"i",&m8, &n8, &l8); // Variable flag

  for(i=0; i<*istk(l1); i++)
    {
      *(dstk(l7)+i) = out_hp[i];
    } 
  *istk(l8) = out_flag;

  free(out_hp);

  LhsVar(1) = 8;
  LhsVar(2) = 9;
  LhsVar(3) = 10;

  PutLhsVar();

  return 0;
}

////////////////////
// interface endp //
////////////////////

int interface_endp_(char * fname)
{
  int l1, m1, n1; // variable n
  int l2, m2, n2; // variable x(n)
  int l3, m3, n3; // variable l(n)
  int l4, m4, n4; // variable u(n)
  int l5, m5, n5; // variable m
  int l6, m6, n6; // variable lambda(m)
  int l7, m7, n7; // variable rho(m)
  int l8, m8, n8; // variable equatn(m)
  int l9, m9, n9; // variable linear(m)
  int l10, m10, n10; // variable function name

  int ierr = 0;

  GetRhsVar(1,"c",&m10, &n10, &l10); // Variable function name
  GetRhsVar(2,"i",&m1, &n1, &l1);    // Variable n
  GetRhsVar(3,"d",&m2, &n2, &l2);    // Variable x(n)
  GetRhsVar(4,"d",&m3, &n3, &l3);    // Variable l(n)
  GetRhsVar(5,"d",&m4, &n4, &l4);    // Variable u(n)
  GetRhsVar(6,"i",&m5, &n5, &l5);    // Variable m
  GetRhsVar(7,"d",&m6, &n6, &l6);    // Variable lambda(m)
  GetRhsVar(8,"d",&m7, &n7, &l7);    // Variable rho(m)
  GetRhsVar(9,"i",&m8, &n8, &l8);    // Variable equatn(m)
  GetRhsVar(10,"i",&m9, &n9, &l9);   // Variable linear(m)

  if (strcmp(cstk(l10),"america")==0)
    {
      ierr = america_endp__(istk(l1),  // Variable n
			    dstk(l2),  // Variable x(n)
			    dstk(l3),  // Variable l(n)
			    dstk(l4),  // Variable u(n)
			    istk(l5),  // Variable m
			    dstk(l6),  // Variable lambda(m)
			    dstk(l7),  // Variable rho(m)
			    lstk(l8),  // Variable equatn(m)
			    lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"bratu3db")==0)
    {
      ierr = bratu3db_endp__(istk(l1),  // Variable n
			     dstk(l2),  // Variable x(n)
			     dstk(l3),  // Variable l(n)
			     dstk(l4),  // Variable u(n)
			     istk(l5),  // Variable m
			     dstk(l6),  // Variable lambda(m)
			     dstk(l7),  // Variable rho(m)
			     lstk(l8),  // Variable equatn(m)
			     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"cache")==0)
    {
      ierr = cache_endp__(istk(l1),  // Variable n
			  dstk(l2),  // Variable x(n)
			  dstk(l3),  // Variable l(n)
			  dstk(l4),  // Variable u(n)
			  istk(l5),  // Variable m
			  dstk(l6),  // Variable lambda(m)
			  dstk(l7),  // Variable rho(m)
			  lstk(l8),  // Variable equatn(m)
			  lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"condor")==0)
    {
      ierr = condor_endp__(istk(l1),  // Variable n
			   dstk(l2),  // Variable x(n)
			   dstk(l3),  // Variable l(n)
			   dstk(l4),  // Variable u(n)
			   istk(l5),  // Variable m
			   dstk(l6),  // Variable lambda(m)
			   dstk(l7),  // Variable rho(m)
			   lstk(l8),  // Variable equatn(m)
			   lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"contor2")==0)
    {
      ierr = contor2_endp__(istk(l1),  // Variable n
			    dstk(l2),  // Variable x(n)
			    dstk(l3),  // Variable l(n)
			    dstk(l4),  // Variable u(n)
			    istk(l5),  // Variable m
			    dstk(l6),  // Variable lambda(m)
			    dstk(l7),  // Variable rho(m)
			    lstk(l8),  // Variable equatn(m)
			    lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"contor")==0)
    {
      ierr = contor_endp__(istk(l1),  // Variable n
			   dstk(l2),  // Variable x(n)
			   dstk(l3),  // Variable l(n)
			   dstk(l4),  // Variable u(n)
			   istk(l5),  // Variable m
			   dstk(l6),  // Variable lambda(m)
			   dstk(l7),  // Variable rho(m)
			   lstk(l8),  // Variable equatn(m)
			   lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"ellipsoid")==0)
    {
      ierr = ellipsoid_endp__(istk(l1),  // Variable n
			      dstk(l2),  // Variable x(n)
			      dstk(l3),  // Variable l(n)
			      dstk(l4),  // Variable u(n)
			      istk(l5),  // Variable m
			      dstk(l6),  // Variable lambda(m)
			      dstk(l7),  // Variable rho(m)
		  	      lstk(l8),  // Variable equatn(m)
			      lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"genpack-cc-mina")==0)
    {
      ierr = genpack_cc_mina_endp__(istk(l1),  // Variable n
				    dstk(l2),  // Variable x(n)
				    dstk(l3),  // Variable l(n)
				    dstk(l4),  // Variable u(n)
				    istk(l5),  // Variable m
				    dstk(l6),  // Variable lambda(m)
				    dstk(l7),  // Variable rho(m)
				    lstk(l8),  // Variable equatn(m)
				    lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"genpack-csq-mina")==0)
    {
      ierr = genpack_csq_mina_endp__(istk(l1),  // Variable n
				     dstk(l2),  // Variable x(n)
				     dstk(l3),  // Variable l(n)
				     dstk(l4),  // Variable u(n)
				     istk(l5),  // Variable m
				     dstk(l6),  // Variable lambda(m)
				     dstk(l7),  // Variable rho(m)
				     lstk(l8),  // Variable equatn(m)
				     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"hardcube")==0)
    {
      ierr = hardcube_endp__(istk(l1),  // Variable n
			     dstk(l2),  // Variable x(n)
			     dstk(l3),  // Variable l(n)
			     dstk(l4),  // Variable u(n)
			     istk(l5),  // Variable m
			     dstk(l6),  // Variable lambda(m)
			     dstk(l7),  // Variable rho(m)
			     lstk(l8),  // Variable equatn(m)
			     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"hardspheres")==0)
    {
      ierr = hardspheres_endp__(istk(l1),  // Variable n
				dstk(l2),  // Variable x(n)
				dstk(l3),  // Variable l(n)
				dstk(l4),  // Variable u(n)
				istk(l5),  // Variable m
				dstk(l6),  // Variable lambda(m)
				dstk(l7),  // Variable rho(m)
				lstk(l8),  // Variable equatn(m)
				lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"kissing2")==0)
    {
      ierr = kissing2_endp__(istk(l1),  // Variable n
			     dstk(l2),  // Variable x(n)
			     dstk(l3),  // Variable l(n)
			     dstk(l4),  // Variable u(n)
			     istk(l5),  // Variable m
			     dstk(l6),  // Variable lambda(m)
			     dstk(l7),  // Variable rho(m)
			     lstk(l8),  // Variable equatn(m)
			     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"kissing")==0)
    {
      ierr = kissing_endp__(istk(l1),  // Variable n
			    dstk(l2),  // Variable x(n)
			    dstk(l3),  // Variable l(n)
			    dstk(l4),  // Variable u(n)
			    istk(l5),  // Variable m
			    dstk(l6),  // Variable lambda(m)
			    dstk(l7),  // Variable rho(m)
			    lstk(l8),  // Variable equatn(m)
			    lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"location")==0)
    {
      ierr = location_endp__(istk(l1),  // Variable n
			     dstk(l2),  // Variable x(n)
			     dstk(l3),  // Variable l(n)
			     dstk(l4),  // Variable u(n)
			     istk(l5),  // Variable m
			     dstk(l6),  // Variable lambda(m)
			     dstk(l7),  // Variable rho(m)
			     lstk(l8),  // Variable equatn(m)
			     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"mountain1")==0)
    {
      ierr = mountain1_endp__(istk(l1),  // Variable n
			      dstk(l2),  // Variable x(n)
			      dstk(l3),  // Variable l(n)
			      dstk(l4),  // Variable u(n)
			      istk(l5),  // Variable m
			      dstk(l6),  // Variable lambda(m)
			      dstk(l7),  // Variable rho(m)
			      lstk(l8),  // Variable equatn(m)
			      lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"mountain2")==0)
    {
      ierr = mountain2_endp__(istk(l1),  // Variable n
			      dstk(l2),  // Variable x(n)
			      dstk(l3),  // Variable l(n)
			      dstk(l4),  // Variable u(n)
			      istk(l5),  // Variable m
			      dstk(l6),  // Variable lambda(m)
			      dstk(l7),  // Variable rho(m)
			      lstk(l8),  // Variable equatn(m)
			      lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"packccmn")==0)
    {
      ierr = packccmn_endp__(istk(l1),  // Variable n
			     dstk(l2),  // Variable x(n)
			     dstk(l3),  // Variable l(n)
			     dstk(l4),  // Variable u(n)
			     istk(l5),  // Variable m
			     dstk(l6),  // Variable lambda(m)
			     dstk(l7),  // Variable rho(m)
			     lstk(l8),  // Variable equatn(m)
			     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"packccmn-feas")==0)
    {
      ierr = packccmn_feas_endp__(istk(l1),  // Variable n
				  dstk(l2),  // Variable x(n)
				  dstk(l3),  // Variable l(n)
				  dstk(l4),  // Variable u(n)
				  istk(l5),  // Variable m
				  dstk(l6),  // Variable lambda(m)
				  dstk(l7),  // Variable rho(m)
				  lstk(l8),  // Variable equatn(m)
				  lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"packcrmn-feas")==0)
    {
      ierr = packcrmn_feas_endp__(istk(l1),  // Variable n
				  dstk(l2),  // Variable x(n)
				  dstk(l3),  // Variable l(n)
				  dstk(l4),  // Variable u(n)
				  istk(l5),  // Variable m
				  dstk(l6),  // Variable lambda(m)
				  dstk(l7),  // Variable rho(m)
				  lstk(l8),  // Variable equatn(m)
				  lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"pedazos4")==0)
    {
      ierr = pedazos4_endp__(istk(l1),  // Variable n
			     dstk(l2),  // Variable x(n)
			     dstk(l3),  // Variable l(n)
			     dstk(l4),  // Variable u(n)
			     istk(l5),  // Variable m
			     dstk(l6),  // Variable lambda(m)
			     dstk(l7),  // Variable rho(m)
			     lstk(l8),  // Variable equatn(m)
			     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"piecefit")==0)
    {
      ierr = piecefit_endp__(istk(l1),  // Variable n
			     dstk(l2),  // Variable x(n)
			     dstk(l3),  // Variable l(n)
			     dstk(l4),  // Variable u(n)
			     istk(l5),  // Variable m
			     dstk(l6),  // Variable lambda(m)
			     dstk(l7),  // Variable rho(m)
			     lstk(l8),  // Variable equatn(m)
			     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"simfock2")==0)
    {
      ierr = simfock2_endp__(istk(l1),  // Variable n
			     dstk(l2),  // Variable x(n)
			     dstk(l3),  // Variable l(n)
			     dstk(l4),  // Variable u(n)
			     istk(l5),  // Variable m
			     dstk(l6),  // Variable lambda(m)
			     dstk(l7),  // Variable rho(m)
			     lstk(l8),  // Variable equatn(m)
			     lstk(l9)); // Variable linear(m)
    }
  else if (strcmp(cstk(l10),"simfock")==0)
    {
      ierr = simfock_endp__(istk(l1),  // Variable n
			    dstk(l2),  // Variable x(n)
			    dstk(l3),  // Variable l(n)
			    dstk(l4),  // Variable u(n)
			    istk(l5),  // Variable m
			    dstk(l6),  // Variable lambda(m)
			    dstk(l7),  // Variable rho(m)
			    lstk(l8),  // Variable equatn(m)
			    lstk(l9)); // Variable linear(m)
    }
  else
    {
      Scierror(999,"Wrong function name\n");
      return 0;
    }

  return 0;
}
