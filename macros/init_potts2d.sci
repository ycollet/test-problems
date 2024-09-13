function Var = init_potts2d(Size, Level)
Var = zeros(Size,Size);
for i=1:Size
  for j=1:Size
    Value = floor(rand(1,1)*Level);
    if (rand(1,1)<0.5) then
      Var(i,j) = -Value;
    else
      Var(i,j) = Value;
    end
  end
end
endfunction
