function x_next = f_cvturn(x, T)
% x = [x;y;v;psi;w]
px=x(1); py=x(2); v=x(3); psi=x(4); w=x(5);

x_next = [ px + v*cos(psi)*T;
           py + v*sin(psi)*T;
           v;
           normalizeAngle(psi + w*T);
           w ];
end