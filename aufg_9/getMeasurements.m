function Z = getMeasurements(x_targets)

var =[2 2 0.05]';

%measurement z = (x,y,psi)
%state: X = (x,y,psi,v, vpsi)

p = randperm(numel(x_targets));
for i=1:numel(p)
  x = x_targets{p(i)};
  M = [1 0 0;
       0 1 0;
       0 0 1];
  H = [M zeros(size(M,1),size(x,1)-size(M,2));];

  z = H*x + randn(size(M,1),1).*sqrt(var);
  Z(i).tid = p(i);
  Z(i).z = z;
end
