function Res = NK(Var, N, K, Seed)
// This function is the NK test function for combinatorial optimization
// N : the size of the binary code
// K : the complexity of the problem ((2*K+1) must be inferior or equal to N). 
//     If K = 0, the problem is simple. If K = (N-1)/2, the problem is complex.
// Seed : a parameter to set the rand seed

if (type(Var)==10) then
  AuxString = Var;
  Var = zeros(1,length(AuxString));
  for i=1:length(AuxString)
    Var(i) = eval(part(AuxString,i));
  end
end
if (size(Var,1)==1) then
  Var = Var';
end
if (size(Var,1)~=N) then
  error('size of Var must be equal to N");
end
if ((2*K+1)>N) then
  error('(2*K+1) must be inferior or equal to N");
end
if (~isdef('Seed','local')) then
  Seed = 0;
end

Aux_Var = [Var; Var; Var];
Res     = 0;

for i=1:N
  Computation = 0;
  for j=-K:K
    Computation = Computation + Aux_Var(size(Var,1)+j)*2^(j+K);
  end
  rand('seed',Seed + Computation);
  Res = Res + rand(1,1);
end
Res = Res / (2*K+1);
endfunction
