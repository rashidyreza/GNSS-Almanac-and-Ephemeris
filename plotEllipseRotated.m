function plotEllipseRotated(a,b,C,angle)

% % % % text(C(1,1)+10,C(1,2)+20,kk,'Color','red','FontSize',10);
% % % % text(C(1,1)+30,C(1,2)+10,Angles,'Color','red','FontSize',10);


% range to plot over
%------------------------------------
N = 100;hold on;
theta = 0:1/N:2*pi+1/N;
% Parametric equation of the ellipse
for i=1:size(angle,1)
    phi = angle(i);
    
    %----------------------------------------
    ellipse_x_r  = a(i)*cos( theta );
    ellipse_y_r  = b(i)*sin( theta );
    
    %Define a rotation matrix
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    
    %let's rotate the ellipse to some angle phi
    r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
    
    % Coordinate transform (since your ellipse is axis aligned)
    %----------------------------------------
    X = r_ellipse';
    X(1,:) = X(1,:) + C(i,1);
    X(2,:) = X(2,:) + C(i,2);
    % Plot     %----------------------------------------
    plot(X(1,:),X(2,:));
    
    plot(C(i,1),C(i,2),'.');
    axis equal;
    grid;
%         for j=i+1:size(angle,1);
%            line([C(i,1);C(j,1)],[C(i,2);C(j,2)]); 
%         end
end
