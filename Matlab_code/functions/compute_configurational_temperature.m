function cft = compute_configurational_temperature(X, dV)

f = -dV(X);

[d, nr] = size(X);
cft=0;
for i=1:nr
    
   cft = cft - X(:,i).' * f(:,i); 
end
cft = cft / nr / d;


end