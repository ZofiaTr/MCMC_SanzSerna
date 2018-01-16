function res  = dV(x, k, d)


% radius
r = @(x) norm(x, 2);
% gradient of the potential

if r(x) ==0
    res=zeros(size(x));
  else  
    res =k * (r(x) - 1).*x./r(x);
end
end