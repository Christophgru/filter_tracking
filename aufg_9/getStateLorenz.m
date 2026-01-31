function x_new = getStateLorenz(x_old,T,tid)
timeStepsMeasure = 3;
dt = T / timeStepsMeasure;
x_new = x_old;
for i=1:timeStepsMeasure
  x_new = getStateLorenzIter(x_new, dt, tid);
end



function x_new = getStateLorenzIter(x_old, dt, tid)

persistent x;
persistent v_vec;

if isempty(x)
  x(tid)=4;
end
if numel(x)<tid
  x(tid)=x(tid-1)+2;
end

if isempty(v_vec) || numel(v_vec)<tid
  v_vec{tid} = [0 0 0]'; 
end


y = x_old(1);
z = x_old(2);


%K = [28,46.92,4];
K = [8,46.92,2];

vx = K(1)*(y-x(tid));
vy = x(tid)*(K(2)-z)-y;   
vz = x(tid)*y-K(3)*z;

x(tid) = x(tid)+dt*vx;

y = y+dt*vy;
z = z+dt*vz;

v = sqrt(vy^2+vz^2);
psi = atan2(vz, vy);
psi_old = atan2(v_vec{tid}(2), v_vec{tid}(1));
vpsi = normalizeAngle(psi-psi_old)/dt;


x_new = [y z psi v vpsi]';

v_vec{tid} = [vx vy vz]';


function phi = normalizeAngle(phi_in)
% correct the angular value difference if one and only one value is < 0:
% the angular measurement range is (-180:180], so the difference between
% the values are different than the numeric result!!
phi = abs(phi_in);
% normalize phi to value in interval [-pi;pi]
phi = mod(phi, 2*pi);
if(phi > pi)
  phi = phi - 2*pi;
end
phi = phi*sign(phi_in);