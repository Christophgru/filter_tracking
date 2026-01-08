
function nis = calcNIS(nu, S)
nu = nu(:);
nis = nu' * (S \ nu);
end
