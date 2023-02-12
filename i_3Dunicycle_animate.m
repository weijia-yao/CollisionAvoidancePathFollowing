%% plot the desired paths
figure; hold on; grid on; set(gcf,'color','w'); axis equal;
[x, y, z]=sphere;
h1=surf(x-2.8,y,z, 'FaceAlpha',0.2, 'EdgeColor', 'none', 'FaceColor', [0.8 0.5 0]);  % reactive area
hold on;
plot3(-2.8,0,0,'r.','MarkerSize',100);    % obstacle

syms x y
e11 = 0.035*log(((x + 3.0)^2 + y^2)^(1/2) + 1.0)*((x + 3.0)^2 + y^2) - 0.048*log(((x - 1.5)^2 + y^2)^(1/2) + 1.0)*((x - 1.5)^2 + y^2) - 0.048*log(((y - 1.3)^2 + (x + 0.75)^2)^(1/2) + 1.0)*((y - 1.3)^2 + (x + 0.75)^2) - 0.048*log(((y + 1.3)^2 + (x + 0.75)^2)^(1/2) + 1.0)*((y + 1.3)^2 + (x + 0.75)^2) + 0.035*log(((y - 2.6)^2 + (x - 1.5)^2)^(1/2) + 1.0)*((y - 2.6)^2 + (x - 1.5)^2) + 0.035*log(((y + 2.6)^2 + (x - 1.5)^2)^(1/2) + 1.0)*((y + 2.6)^2 + (x - 1.5)^2) - 1.0;
fimplicit(e11, 'Color', 'red', 'LineWidth', 2) % desired path
% legend([h1, h2],'reactive area', 'desired path');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(-140,20)
hAxes = gca;

%% plot starting position
plot3(p(1,1),p(1,2),p(1,3),'bo','LineWidth', 1);

h = animatedline(hAxes);
hh = animatedline(hAxes, 'Color','blue', 'Marker', '.', 'MarkerSize', 20);
htri = animatedline(hAxes, 'MaximumNumPoints', Inf);

%% animation
numpoints = size(p,1);
drawnum = 10;
for i = 1 : drawnum : numpoints-(drawnum-1)
    % dis = sprintf("Time (s): %0.3f", e.Time(i+(drawnum-1)));
    dis = sprintf("Time (s): %0.3f", tout(i+(drawnum-1)));
    tt = text(hAxes,1,1,2, dis);
    
    % trajectory 
    addpoints(h, p(i:i+(drawnum-1),1), p(i:i+(drawnum-1),2), p(i:i+(drawnum-1),3));
    % trajectory point
    addpoints(hh, p(i+(drawnum-1),1), p(i+(drawnum-1),2), p(i+(drawnum-1),3));
    
    if(IS_UNICYCLE)
        theta = p(i+(drawnum-1),4);
        p1 = [p(i+(drawnum-1),1); p(i+(drawnum-1),2); p(i+(drawnum-1),3)];
        vel = [dpdt(i+(drawnum-1),1); dpdt(i+(drawnum-1),2); dpdt(i+(drawnum-1),3)];      
        L = 0.5; 
        p1 = p1 - 0.433*L*vel/norm(vel);
        tript = cal_tri_pts(p1, vel, L);		% get four points to draw a triangle and then draw it
        addpoints(htri, tript(1,3), tript(2,3), tript(3,3));
        addpoints(htri, tript(1,1), tript(2,1), tript(3,1));
        addpoints(htri, tript(1,2), tript(2,2), tript(3,2));
        addpoints(htri, tript(1,3), tript(2,3), tript(3,3));
        addpoints(htri, tript(1,4), tript(2,4), tript(3,4));
        addpoints(htri, tript(1,1), tript(2,1), tript(3,1));
    end
    
    drawnow
    pause(0.01) 
    
    delete(tt)
    clearpoints(hh)
    if(IS_UNICYCLE)
        clearpoints(htri)
    end
end

function p = cal_tri_pts(p1, v, L)
% input:    p1 -- the midpoint of one edge of a triangle; COLUMN vector
%           v  -- a velocity vector; pointing from p1 to a vertex
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