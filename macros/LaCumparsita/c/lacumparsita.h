//
// Functions *_inip__
//

typedef double doublereal;
//typedef int integer;
typedef int logical;

extern int america_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int bratu3db_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *nppar, doublereal *seed);
extern int cache_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int condor_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *n_in__, integer *p_in__, integer *q_in__);
extern int contor2_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *n_in__);
extern int contor_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *n_in__);
extern int ellipsoid_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *nd_in__, integer *np_in__);
extern int genpack_cc_mina_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, doublereal *iterad_in__, integer *nite_in__, doublereal *seed_in__);
extern int genpack_csq_mina_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, doublereal *iterad_in__, integer *nite_in__, doublereal *seed_in__);
extern int hardcube_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *nd_in__, integer *np_in__, doublereal *seed_in__);
extern int hardspheres_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *nd_in__, integer *np_in__, doublereal *seed_in__);
extern int kissing2_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *nd_in__, integer *np_in__, doublereal *seed_in__);
extern int kissing_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *nd_in__, integer *np_in__, doublereal *seed_in__);
extern int location_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int mountain1_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *np_in__, doublereal *dmax2_in__, doublereal *seed_in__);
extern int mountain2_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *np_in__, doublereal *dmax2_in__, doublereal *seed_in__);
extern int packccmn_feas_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *nd_in__, integer *nite_in__, doublereal *iterad_in__, doublereal *objrad_in__, doublereal *seed_in__);
extern int packccmn_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int packcrmn_feas_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *nd_in__, integer *nite_in__, doublereal *iterad_in__, doublereal *objdim_in__, doublereal *seed_in__);
extern int pedazos4_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int piecefit_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *n_in__, integer *p_in__);
extern int simfock2_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *n_in__, integer *k_in__);
extern int simfock_inip__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear, integer *n_in__, integer *k_in__);

//
// Functions *_evalf__
//
extern int america_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int bratu3db_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int cache_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int condor_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int contor2_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int contor_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int ellipsoid_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int genpack_cc_mina_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int genpack_csq_mina_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int hardcube_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int hardspheres_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int kissing2_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int kissing_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int location_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int mountain1_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int mountain2_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int packccmn_feas_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int packccmn_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int packcrmn_feas_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int pedazos4_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int piecefit_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int simfock2_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);
extern int simfock_evalf__(integer *n, doublereal *x, doublereal *f, integer *flag__);

//
// Functions *_evalg__
//
extern int america_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int bratu3db_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int cache_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int condor_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int contor2_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int contor_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int ellipsoid_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int genpack_cc_mina_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int genpack_csq_mina_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int hardcube_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int hardspheres_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int kissing2_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int kissing_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int location_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int mountain1_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int mountain2_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int packccmn_feas_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int packccmn_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int packcrmn_feas_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int pedazos4_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int piecefit_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int simfock2_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);
extern int simfock_evalg__(integer *n, doublereal *x, doublereal *g, integer *flag__);

//
// Functions *_evalh__
//
extern int america_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int bratu3db_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *hnnz, integer *flag__);
extern int cache_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int condor_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int contor2_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int contor_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int ellipsoid_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *hnnz, integer *flag__);
extern int genpack_cc_mina_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int genpack_csq_mina_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int hardcube_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int hardspheres_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *hnnz, integer *flag__);
extern int kissing2_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *hnnz, integer *flag__);
extern int kissing_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *hnnz, integer *flag__);
extern int location_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *hnnz, integer *flag__);
extern int mountain1_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int mountain2_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int packccmn_feas_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *hnnz, integer *flag__);
extern int packccmn_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int packcrmn_feas_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *hnnz, integer *flag__);
extern int pedazos4_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int piecefit_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int simfock2_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);
extern int simfock_evalh__(integer *n, doublereal *x, integer *hlin, integer *hcol, doublereal *hval, integer *nnzh, integer *flag__);

//
// Functions *_evalc__
//
extern int america_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int bratu3db_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int cache_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int condor_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int contor2_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int contor_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int ellipsoid_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int genpack_cc_mina_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int genpack_csq_mina_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int hardcube_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int hardspheres_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int kissing2_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int kissing_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int location_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int mountain1_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int mountain2_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int packccmn_feas_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int packccmn_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int packcrmn_feas_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int pedazos4_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int piecefit_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int simfock2_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);
extern int simfock_evalc__(integer *n, doublereal *x, integer *ind, doublereal *c__, integer *flag__);

//
// Functions *_evaljac__
//
extern int america_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int bratu3db_evaljac__(integer *n, doublereal *x, integer *ind, integer *jcvar, doublereal *jcval, integer *jcnnz, integer *flag__);
extern int cache_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int condor_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int contor2_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int contor_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int ellipsoid_evaljac__(integer *n, doublereal *x, integer *ind, integer *jcvar, doublereal *jcval, integer *jcnnz, integer *flag__);
extern int genpack_cc_mina_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int genpack_csq_mina_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int hardcube_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int hardspheres_evaljac__(integer *n, doublereal *x, integer *ind, integer *jcvar, doublereal *jcval, integer *jcnnz, integer *flag__);
extern int kissing2_evaljac__(integer *n, doublereal *x, integer *ind, integer *jcvar, doublereal *jcval, integer *jcnnz, integer *flag__);
extern int kissing_evaljac__(integer *n, doublereal *x, integer *ind, integer *jcvar, doublereal *jcval, integer *jcnnz, integer *flag__);
extern int location_evaljac__(integer *n, doublereal *x, integer *ind, integer *jcvar, doublereal *jcval, integer *jcnnz, integer *flag__);
extern int mountain1_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int mountain2_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int packccmn_feas_evaljac__(integer *n, doublereal *x, integer *ind, integer *jcvar, doublereal *jcval, integer *jcnnz, integer *flag__);
extern int packccmn_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int packcrmn_feas_evaljac__(integer *n, doublereal *x, integer *ind, integer *jcvar, doublereal *jcval, integer *jcnnz, integer *flag__);
extern int pedazos4_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int piecefit_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int simfock2_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);
extern int simfock_evaljac__(integer *n, doublereal *x, integer *ind, integer *indjac, doublereal *valjac, integer *nnzjac, integer *flag__);

//
// Functions *_evalhc__
//
extern int america_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int bratu3db_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int cache_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int condor_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int contor2_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int contor_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int ellipsoid_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *hcnnz, integer *flag__);
extern int genpack_cc_mina_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int genpack_csq_mina_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int hardcube_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int hardspheres_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *hcnnz, integer *flag__);
extern int kissing2_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *hcnnz, integer *flag__);
extern int kissing_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *hcnnz, integer *flag__);
extern int location_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *hcnnz, integer *flag__);
extern int mountain1_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int mountain2_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int packccmn_feas_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *hcnnz, integer *flag__);
extern int packccmn_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int packcrmn_feas_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *hcnnz, integer *flag__);
extern int pedazos4_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int piecefit_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int simfock2_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);
extern int simfock_evalhc__(integer *n, doublereal *x, integer *ind, integer *hclin, integer *hccol, doublereal *hcval, integer *nnzhc, integer *flag__);

//
// Functions *_evalhlp__
//
extern int america_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int bratu3db_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int cache_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int condor_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int contor2_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int contor_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int ellipsoid_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int genpack_cc_mina_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int genpack_csq_mina_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int hardcube_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int hardspheres_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int kissing2_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int kissing_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int location_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int mountain1_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int mountain2_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int packccmn_feas_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int packccmn_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int packcrmn_feas_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int pedazos4_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int piecefit_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int simfock2_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);
extern int simfock_evalhlp__(integer *n, doublereal *x, integer *m, doublereal *lambda, doublereal *p, doublereal *hp, logical *goth, integer *flag__);

//
// Functions *_endp__
//
extern int america_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int bratu3db_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int cache_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int condor_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int contor2_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int contor_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int ellipsoid_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int genpack_cc_mina_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int genpack_csq_mina_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int hardcube_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int hardspheres_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int kissing2_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int kissing_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int location_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int mountain1_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int mountain2_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int packccmn_feas_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int packccmn_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int packcrmn_feas_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int pedazos4_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int piecefit_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int simfock2_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
extern int simfock_endp__(integer *n, doublereal *x, doublereal *l, doublereal *u, integer *m, doublereal *lambda, doublereal *rho, logical *equatn, logical *linear);
