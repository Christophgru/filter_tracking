function x5 = trueToFilterState(x_true)
% Map simulator true state to EKF state: [x;y;v;psi;w]
x_true = x_true(:);

x5 = zeros(5,1);
x5(1) = x_true(1);
x5(2) = x_true(2);

% you use x_true(3) as speed in plotting -> keep that
x5(3) = x_true(3);

% you use x_true(4) as heading in plotting -> keep that
x5(4) = x_true(4);

% if simulator provides w as 5th component, use it, else 0
if length(x_true) >= 5
    x5(5) = x_true(5);
else
    x5(5) = 0;
end
end
