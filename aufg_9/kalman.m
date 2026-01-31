function kf = kalman(action, p1, p2, p3, p4, p5)
%  kf(nrKF) = kalman('init', z1, z2, T);
%  kf(nrKF) = kalman('predict', kf(nrKF), T, qv, qw); 
%  kf(nrKF) = kalman('update', kf(nrKF), Z(zid).z, R, z_pred, H);

switch(action)
  case 'init'
    z1 = p1;
    z2 = p2;
    T = p3;
    R = p4;

    x_est(1:3,1) = z2;
    diff = z2-z1;
    posDiff = diff(1:2);    
    v = sqrt(sum(posDiff.^2))/T;
    w = diff(3) / T;
    x_est(4,1) = v;
    x_est(5,1) = w;
    
    P_est = [R(1,1) 0 0 0 0 ;
      0 R(1,1) 0 0 0 ;
      0 0 R(3,3) 0 0 ;
      0 0 0 2*R(1,1)/T^2 0  ;
      0 0 0 0 2*R(3,3)/T^2];

    
    kf.x_est = x_est;
    kf.P_est = P_est;
    kf.x_pred = [];
    kf.P_pred = [];
    kf.inno.S = [];
    kf.inno.v = [];
    kf.HPH = [];
    kf.K = [];
    kf.z = [];
    kf.z_pred = [];

  case 'predict'
    kf = p1;
    T = p2;
    qv = p3;
    qw = p4;

    x_est = kf.x_est;
    P_est = kf.P_est;
    x_old = x_est(1);
    y_old =  x_est(2);
    phi_old =  x_est(3);
    v_old =  x_est(4);
    w_old =  x_est(5);
    
    
    %Modell: CV/CT
    
    % the process model matrix (non-linear)
    % CV/CT-model: speed and orientation do not change
    % we use the non-linear process model in the
    % prediction of the state estimate
    
    x_pred = [];
    
    x_pred(1,1) = x_old + v_old*T*cos(phi_old);
    x_pred(2,1) = y_old + v_old*T*sin(phi_old);
    x_pred(3,1) = phi_old + T * w_old;
    x_pred(4,1) = v_old;
    x_pred(5,1) = w_old;
    
    x_pred(3,1) = normalizeAngle(x_pred(3,1));
    
    if(w_old>eps)
      % x y phi v w
      % jacobian matrix of F
      Fx = [1 0 -v_old*T*sin(phi_old)  T*cos(phi_old) 0 ;
        0 1 v_old*T*cos(phi_old)  T*sin(phi_old) 0 ;
        0 0 1 0 T ;
        0 0 0 1 0 ;
        0 0 0 0 1;
        ];
    else
      Fx = [1 0 -v_old*T*sin(phi_old)  T*cos(phi_old) 0 ;
        0 1 v_old*T*cos(phi_old)  T*sin(phi_old) 0 ;
        0 0 1 0 0 ;
        0 0 0 1 0 ;
        0 0 0 0 0;
        ];
    end
    
    
    
    % the process noise (cf. skript)
    Qxx = 1/3*cos(phi_old)^2*qv*T^3;
    Qxy = 1/3*cos(phi_old)*qv*sin(phi_old)*T^3;
    Qxv = 1/2 *cos(phi_old)*qv*T^2;
    
    Qyy = 1/3*(sin(phi_old)^2*qv)*T^3;
    Qyv = 1/2 * sin(phi_old)*qv*T^2;
    
    Qvv = qv * T;
    Qpp = 1/3*qw*T^3;
    Qpw = 1/2*T^2*qw;
    Qww = T*qw;
    
    
    
    %     x    y  p    v  w
    Q = [ Qxx Qxy 0   Qxv   0 ;
      Qxy Qyy 0   Qyv   0 ;
      0   0 Qpp   0 Qpw ;
      Qxv Qyv   0 Qvv   0 ;
      0   0 Qpw   0 Qww ;
      ];
    
    
    % the kalman equations....
    P_pred = Fx*P_est*Fx' + Q;
    
    kf.x_pred = x_pred;
    kf.P_pred = P_pred;
    
  case 'update'
  %  kf(nrKF) = kalman('update', kf(nrKF), Z(zid).z, R, z_pred, H);
  kf = p1;
  z = p2;
  R = p3;
  z_pred = p4;
  H = p5;
  
  
  P_pred = kf.P_pred;
  x_pred = kf.x_pred;
  
  HPH = H*P_pred*H';
  S = HPH+R;
  v = z-z_pred;
  v(3) = normalizeAngle(v(3));
  
 
  
  K = P_pred*H'*inv(S);
  kf.x_est = x_pred + K*v;
  kf.P_est = P_pred - K*S*K';
  kf.inno.S = S;
  kf.inno.v = v;
  kf.HPH = HPH;
  kf.K = K;
  kf.z = z;
  kf.z_pred = z_pred;
end