%% plot the desired path
figure; hold on; grid on; set(gcf,'color','w'); axis equal;
[x, y, z]=sphere;
h1=surf(x-2.8,y,z, 'FaceAlpha',0.2, 'EdgeColor', 'none', 'FaceColor', [0.8 0.5 0]);  % reactive area
radius = sqrt(0.28);
h2=surf(x*radius-2.8,y*radius,z*radius, 'FaceAlpha',0.2, 'EdgeColor', 'none', 'FaceColor', [0.8 0 0.5]);  % repulsive area
hold on;
h3=plot3(-2.8,0,0,'r.','MarkerSize',100);    % obstacle

syms x y
e11 = 0.035*log(((x + 3.0)^2 + y^2)^(1/2) + 1.0)*((x + 3.0)^2 + y^2) - 0.048*log(((x - 1.5)^2 + y^2)^(1/2) + 1.0)*((x - 1.5)^2 + y^2) - 0.048*log(((y - 1.3)^2 + (x + 0.75)^2)^(1/2) + 1.0)*((y - 1.3)^2 + (x + 0.75)^2) - 0.048*log(((y + 1.3)^2 + (x + 0.75)^2)^(1/2) + 1.0)*((y + 1.3)^2 + (x + 0.75)^2) + 0.035*log(((y - 2.6)^2 + (x - 1.5)^2)^(1/2) + 1.0)*((y - 2.6)^2 + (x - 1.5)^2) + 0.035*log(((y + 2.6)^2 + (x - 1.5)^2)^(1/2) + 1.0)*((y + 2.6)^2 + (x - 1.5)^2) - 1.0;
h4=fimplicit(e11, 'Color', 'red', 'LineWidth', 3, 'LineStyle', '--') % desired path

%% plot the traj
h5=plot3(p(:,1), p(:,2), p(:,3), 'm', 'LineWidth',2);
plot3(p(1,1), p(1,2), p(1,3), 'bo','MarkerSize', 8);

if(IS_UNICYCLE)
    index = [1, 82, 208, 400, 570, 680, 756, 868];
    L = 0.3;
    for i=1:length(index)
        id = index(i);
        theta = p(id,4);
        point=[p(id,1); p(id,2); p(id,3)];
        vel = [dpdt(id,1); dpdt(id,2); dpdt(id,3)];
        tript = cal_tri_pts(point, vel, L);
        line_pts = [tript(:,3), tript(:,1), tript(:,2), tript(:,3), tript(:,4), tript(:,1)];
        line(line_pts(1,:), line_pts(2,:), line_pts(3,:), 'LineWidth', 2);
    end
end


xlabel('X'); ylabel('Y'); zlabel('Z');
legend([h1, h2, h4, h5],'reactive surface','repulsive surface','desired path', 'trajectory');
view(39,23)
hold off;


%% function to plot triangles
function p = cal_tri_pts(p1, v, L)
% input:    p1 -- the midpoint of one edge of a triangle; COLUMN vector
%           v  -- a velocity vector; pointing from p1 to a vertex; COLUMN vector
%           L  -- edge length of the triangle
% output:   p = [p1, p2, p3, p4] contaning the verteces of the triangle
%           p2, p3 and p4 (COLUMN vectors); p1 is the same as the input;
%
%                          p3
%                          *                     
%                         /|\
%                        / | \
%                       /  |  \
%                      /___|___\ 
%                     p4   p1   p2
%
    p3=p1+0.866*L*v/norm(v);
    solx = v(2)*(1/(v(1)^2 + v(2)^2))^(1/2);
    soly = -v(1)*(1/(v(1)^2 + v(2)^2))^(1/2);
    vv=double([solx;soly;0]);
    p2=p1+0.5*L*vv;
    p4=p1-0.5*L*vv;
    p=[p1, p2, p3, p4];
end