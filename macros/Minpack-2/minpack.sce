//listoffiles = ['daerfj','datrfj','dchqfj','dcpffj','dcprfj','dctsfj','dedffj','deptfg','depths','deptsp','dfdcfj','dfdcjs','dfdcsp','dficfj','dficjs','dficsp', ...
//               'dgdffj','dgl1fg','dgl1hs','dgl1sp','dgl2co','dgl2fg','dgl2hs','dgl2sp','dhhdfj','diacfj','diadfj','diaofj','diarfj','dierfj','dierjs','diersp', ...
//               'dljcfg','dmsabc','dmsafg','dmsahs','dmsasp','dodcfg','dodchs','dodcps','dodcsp','dpjbds','dpjbfg','dpjbhs','dpjbsp','dsfdfj','dsfdjs','dsfdsp', ...
//               'dsfifj','dsfijs','dsfisp','dsscfg','dsschs','dsscsp','saerfj','satrfj','schqfj','scpffj','scprfj','sctsfj','sedffj','septfg','sepths','septsp', ...
//               'sfdcfj','sfdcjs','sfdcsp','sficfj','sficjs','sficsp','sgdffj','sgl1fg','sgl1hs','sgl1sp','sgl2co','sgl2fg','sgl2hs','sgl2sp','shhdfj','siacfj', ...
//               'siadfj','siaofj','siarfj','sierfj','sierjs','siersp','sljcfg','smsabc','smsafg','smsahs','smsasp','sodcfg','sodchs','sodcps','sodcsp','spjbds', ...
//               'spjbfg','spjbhs','spjbsp','ssfdfj','ssfdjs','ssfdsp','ssfifj','ssfijs','ssfisp','ssscfg','ssschs','ssscsp'];

listoffiles = ['daerfj','datrfj','dchqfj','dcpffj','dcprfj','dctsfj','dedffj','deptfg','depths','deptsp','dfdcfj','dfdcjs','dfdcsp','dficfj','dficjs','dficsp', ...
               'dgdffj','dgl1fg','dgl1hs','dgl1sp','dgl2co','dgl2fg','dgl2hs','dgl2sp','dhhdfj','diacfj','diadfj','diaofj','diarfj','dierfj','dierjs','diersp', ...
               'dljcfg','dminfg','dminhs','dminsp','dminxb','dmsabc','dmsafg','dmsahs','dmsasp','dodcfg','dodchs','dodcps','dodcsp','dpjbds','dpjbfg','dpjbhs', ...
               'dpjbsp','dsfdfj','dsfdjs','dsfdsp','dsfifj','dsfijs','dsfisp','dsscfg','dsschs','dsscsp','dsysfv','dsysjs']

if ~isdef('daerfj') then
  id = link('./libminpack-2.so',listoffiles,'f');
end

//	letters 6 and 7   action
// ---------------   ------
//       fj          function, Jacobian, initial approximation
//       js          Jacobian times vector
//       fg          function, gradient, initial approximation
//       hs          Hessian times vector
//       sp          Jacobian sparsity structure
//       bc          Boundary data
//
// MINPACK-2 test problem collection.
//
// dficfj.f dficjs.f dficsp.f --- flow in a channel
// dsfdfj.f dsfdjs.f dsfdsp.f --- swirling flow between disks
// dierfj.f dierjs.f diersp.f --- incompressible elastic rods
// dsfifj.f dsfijs.f dsfisp.f --- solid fule ignition
// dfdcfj.f dfdcjs.f dfdcsp.f --- flow in a driven cavity
// dhhdfj.f                   --- human heart dipole
// dcpffj.f                   --- combustion of propane: full formulation
// dcprfj.f                   --- combustion of propane: reduced formulation
//
// diacfj.f                   --- Isomerization of alpha-pinene: collocation
// diadfj.f                   --- Isomerization of alpha-pinene: direct
// diaofj.f                   --- Isomerization of alpha-pinene: constraints
// diarfj.f                   --- Isomerization of alpha-pinene: residuals
// dctsfj.f                   --- Coating thickness standardization
// dedffj.f                   --- Exponential data fitting
// dgdffj.f                   --- Gaussian data fitting
// datrfj.f                   --- Analysis of thermistor resistance
// daerfj.f                   --- Analysis on an enzyme reaction
// dchqfj.f                   --- Chebychev quadrature
//
// deptfg.f depths.f deptsp.f --- elastic-plastic torsion
// dpjbfg.f dpjbhs.f dpjbsp.f --- pressure distribution in a journal bearing
//          dpjbds.f
// dmsafg.f dmsahs.f dmsasp.f --- minimal surfaces
// dmsabc.f
// dodcfg.f dodchs.f dodcsp.f --- optimal design with composites
// dodcps.f
// dljcfg.f                   --- Lennard-Jones clusters
// dgl1fg.f dgl1hs.f dgl1sp.f --- 1-d Ginzburg-Landau
// dsscfg.f dsschs.f dsscsp.f --- steady-state combustion
// dgl2fg.f dgl2hs.f dgl2sp.f --- 2-d Ginzburg-Landau
//          dgl2co.f
