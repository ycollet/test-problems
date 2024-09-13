function [f,g,h] = lacumparsita(pbname, x)
listofpb = ['america','bratu3Db','cache','condor','contor2','contor', ...
            'ellipsoid','genpack-cc-mina','genpack-csq-mina','hardcube', ...
            'hardspheres','kissing2','kissing','location','montain1', ...
            'montain2','packccmn','packccmn-feas','packcrmn-feas','pedazos4', ...
            'piecefit','simfock2','simfock'];


endfunction

function [n, x, l, u, m, lambda, rho, equatn, linear] = initp()
[n,m] = call('initp','out',[1,1],1,'i',[1,1],5,'i');
[x, l, u, lambda, rho, equatn, linear] = call('initp','out',[n,1],2,'d', ... // x
                                                            [n,1],3,'d', ... // l
                                                            [n,1],4,'d', ... // u
                                                            [m,1],6,'d', ... // lambda
                                                            [m,1],7,'d', ... // rho
                                                            [m,1],8,'d', ... // equatn
                                                            [m,1],9,'i'); // linear
endfunction

function [f,flag] = evalf(n, x);
[f, flag] = call('evalf',1,n,'i',2,x,'d','out',[1,1],3,'d',[1,1],4,'i');
endfunction

function [g,flag] = evalg(n, x);
[g, flag] = call('evalg',1,n,'i',2,x,'d','out',[n,1],3,'d',[1,1],4,'i');
endfunction

function [hlin, hcol, hval ,flag] = evalh(n, x);
[nnzh, flag] = call('evalh',1,n,'i',2,x,'d','out',[1,1],6,'i',[1,1],7,'i');
[hlin, hcol, hval] = call('evalh',1,n,'i',2,x,'d','out',[nnzh,1],3,'i', ...
                                                        [nnzh,1],4,'i', ...
                                                        [nnzh,1],5,'d');
endfunction

function [c,flag] = evalc(n,x,ind)
[c, flag] = call('evalc',1,n,'i',2,x,'d',3,ind,'i','out',[1,1],4,'d',[1,1],5,'i');
endfunction

function [indjac,valjac,flag] = evaljac(n,x,ind)
[nnzjac, flag] = call('evaljac',1,n,'i',2,x,'d',3,ind,'i','out',[1,1],6,'i',[1,1],7,'i');
[indjac,valjac] = call('evaljac',1,n,'i',2,x,'d',3,ind,'i','out',[nnzjac,1],5,'i',[nnzjac,1],6,'d');
endfunction

function [hclin,hccol,hcjac,flag] = evalhc(n,x,ind)
[nnzhc, flag] = call('evalhc',1,n,'i',2,x,'d',3,ind,'i','out',[1,1],7,'i',[1,1],8,'i');
[hclin,hccol,hcval] = call('evalhc',1,n,'i',2,x,'d',3,ind,'i','out',[nnzhc,1],4,'i',[nnzhc,1],5,'i',[nnzhc,1],6,'d');
endfunction

function [hclin,hccol,hcjac,flag] = evalhlp(n,x,m,lambda,p,goth)
[nnzhc, flag] = call('evalhlp',1,n,'i',2,x,'d',3,ind,'i','out',[1,1],7,'i',[1,1],8,'i');
[hp,goth,flag] = call('evalhlp',1,n,'i', ...
                                2,x,'d', ...
                                3,m,'i', ...
                                4,lambda,'d', ...
                                5,p,'d', ...
                                7,goth,'i', ...
                                'out',[n,1],6,'i', ...
                                      [1,1],7,'i', ...
                                      [1,1],8,'i');
endfunction

function [n, x, l, u, m, lambda, rho, equatn, linear] = endp()
[n,m] = call('endp','out',[1,1],1,'i',[1,1],5,'i');
[x, l, u, lambda, rho, equatn, linear] = call('endp',n,1,'i', ... // n
                                                     x,2,'d', ... // x
                                                     l,3,'d', ... // l
                                                     u,4,'d', ... // u
                                                     lambda,6,'d', ... // lambda
                                                     rho,7,'d', ... // rho
                                                     equatn,8,'d', ... // equatn
                                                     linear,9,'i'); // linear
endfunction
