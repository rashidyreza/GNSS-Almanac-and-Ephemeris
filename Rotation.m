function [R]=Rotation(i,j)
%Rotation matrix
if i==3
    R = [cos(j) sin(j) 0 ; -sin(j) cos(j) 0 ; 0 0 1];
else if i==2
        R = [cos(j) 0 -sin(j) ; 0 1 0 ; sin(j) 0 cos(j)];
    else if i==1
            R = [1 0 0 ; 0 cos(j) sin(j) ; 0 -sin(j) cos(j)];
        end
    end
end
