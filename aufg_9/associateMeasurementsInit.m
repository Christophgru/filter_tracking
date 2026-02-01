function A=associateMeasurementsInit(Z1,Z2)
    % Z1(i).z measurement of one object first time step
    % Z2(i).z measurement of one object second time step
    
    
    A = zeros(numel(Z1), numel(Z2));
    % TODO: implement
    for i =1:numel(Z1)
        closest_neighbor_of_i=0;
        closest_neighbor_dist=inf;
    %calculate the closest neighbour j and then set association Matrix A (i,j)
        for j=1:numel(Z2)
            dist=calc_distance(Z1(i).z,Z2(j).z);
            if (dist<closest_neighbor_dist)
                closest_neighbor_of_i=j;
                closest_neighbor_dist=dist;
            end
        end
        A(i,closest_neighbor_of_i)=1;
    end
end

function dist=calc_distance(z1,z2)
    dist = (z1(1) - z2(1))^2 + (z1(2) - z2(2))^2;
end
