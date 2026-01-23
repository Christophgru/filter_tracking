function x_new = getState(x_old,dt)

    v_cv = 10;
    a_ca = 50;
    lower_time = 40;
    upper_time = 80;
    
    persistent counter;
    persistent toggle_acc; % 0=+a_ca, 1=-a_ca

    if isempty(toggle_acc) || isempty(x_old)
        % reset persistent variable
        toggle_acc = 0; 
    end
    if isempty(counter) || isempty(x_old)
        % reset persistent variable
        counter=0; 
    end

    if counter < lower_time
        % object moves with constant velocity 
%         fprintf(1,'cv - counter:%d\n',counter)
        if isempty(x_old)
            y =  0;
            z = 10;
            vy = v_cv;
            vz = 0;
            ay = 0;
            az = 0;
        else
            y =  x_old(1)+x_old(3)*dt;
            z = 10;
            vy = x_old(3);
            vz = 0;
            ay = 0;
            az = 0;
        end
        
        counter = counter + 1;
    else
        % object moves with constant acceleration
        
        counter = counter + 1;
        if toggle_acc == 0
            % accelerate
            if isempty(x_old)
                y =  0;
                z = 10;
                vy = a_ca*dt;
                vz = 0;
                ay = a_ca;
                az = 0;
            else
                y =  x_old(1)+x_old(3)*dt+0.5*dt^2*a_ca;
                z = 10;
                vy = a_ca*dt+x_old(3);
                vz = 0;
                ay = a_ca;
                az = 0;
            end
            
            if counter == upper_time
            	toggle_acc = 1; 
                counter = 0;
            end
        else
            % brake
            if isempty(x_old)
                y =  0;
                z = 10;
                vy = -a_ca*dt;
                vz = 0;
                ay = -a_ca;
                az = 0;
            else
                y =  x_old(1)+x_old(3)*dt-0.5*dt^2*a_ca;
                z = 10;
                vy = -a_ca*dt+x_old(3);
                vz = 0;
                ay = -a_ca;
                az = 0;
            end
            if counter == upper_time
               toggle_acc = 0; 
               counter = 0;
            end
        end
    end
    x_new = [y z vy vz ay az]';

end

