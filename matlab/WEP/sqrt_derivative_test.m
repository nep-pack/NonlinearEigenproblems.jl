function val = sqrt_derivative_test(a,b,c, d_vec, x)

syms z
f(z) = sqrt( a*z^2 + b*z + c);

d_max = max(d_vec);

val = NaN(length(d_vec),1);
val_idx = 0;

for d = 0:d_max
  idx = find(d_vec==d, 1);
  if(~isempty(idx))
    val_idx = val_idx + 1;
    val(val_idx) = vpa(subs(f,x), 20);
  end
  f = diff(f);
end

clear syms
end
