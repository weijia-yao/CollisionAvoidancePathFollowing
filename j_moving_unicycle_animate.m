WRITE_VIDEO = 0;
WINDOWS = 0;
fontsize = 20;

%% plot the desired paths
figure('Position', [100 100 1024 768]); axis equal; set(gcf,'color','w'); hold on; axis equal; % axis off
syms x y
phi(x,y) = y - sin(x);   % desired path
fimplicit(phi, '--', 'Color', 'red', 'LineWidth', 2)

plot(p(1,1),p(1,2),'bo','LineWidth', 1); %%% plot starting position
xlabel('X'); ylabel('Y');
set(findall(gcf,'-property','FontSize'),'FontSize', fontsize);  % set font size
hAxes = gca; 


%% animation
h = animatedline(hAxes, 'Color', 'magenta', 'LineWidth', 2);
hh = animatedline(hAxes, 'Color','blue', 'Marker', '.', 'MarkerSize', 20, 'LineWidth', 2);
htri = animatedline(hAxes, 'MaximumNumPoints',Inf);
axis([-8 3 -3 3]);

numpoints = size(p,1);
drawnum = 1;
if(WRITE_VIDEO)
    F(numpoints) = struct('cdata',[],'colormap',[]);
    if(WINDOWS) % used in Windows or Mac 
        v = VideoWriter('switch_collision', 'MPEG-4');  
    else
        v = VideoWriter('switch_collision.avi');
    end
    open(v);
end

try
    for i = 1 : drawnum : numpoints-(drawnum-1)
        t = tout(i+(drawnum-1));
        dis = sprintf("Time (s): %0.3f", t);
        tt = text(hAxes, -4, 3, dis);

        addpoints(h, p(i:i+(drawnum-1),1), p(i:i+(drawnum-1),2));  % trajectory 
        addpoints(hh, p(i+(drawnum-1),1), p(i+(drawnum-1),2)); % trajectory point
        
        f = (x+5-obs_v*t)^2/a^2 + y^2/b^2 - 1;  % reactive boundary 
        reacb = fimplicit(hAxes, f, 'Color', 'green', 'LineWidth', 2);
        f = (x+5-obs_v*t)^2/a^2 + y^2/b^2 - 1 - c;  % repulsive boundary 
        repulb = fimplicit(hAxes, f, '-.', 'Color', 'black', 'LineWidth', 2);

        if(IS_UNICYCLE)
            theta = p(i+(drawnum-1),3);
            p1 = [p(i+(drawnum-1),1); p(i+(drawnum-1),2); 0];
            vel = [speed*cos(theta); speed*sin(theta); 0];
            L = 0.2; 
            p1 = p1 - 0.433*L*vel/norm(vel);
            tript = cal_tri_pts(p1, vel, L);		% get four points to draw a triangle and then draw it
            addpoints(htri, tript(1,3), tript(2,3), tript(3,3));
            addpoints(htri, tript(1,1), tript(2,1), tript(3,1));
            addpoints(htri, tript(1,2), tript(2,2), tript(3,2));
            addpoints(htri, tript(1,3), tript(2,3), tript(3,3));
            addpoints(htri, tript(1,4), tript(2,4), tript(3,4));
            addpoints(htri, tript(1,1), tript(2,1), tript(3,1));
        end
        
        if(WRITE_VIDEO)
            cdata = print('-f1', '-RGBImage','-r0');
            F(drawnum) = im2frame(cdata);
            writeVideo(v, F(drawnum));
        else
            drawnow;
            pause(0.01) ;
        end

        delete(tt); delete(reacb); delete(repulb);
        clearpoints(hh)
        if(IS_UNICYCLE)
            clearpoints(htri)
        end
    end
catch draw_error
    disp('Drawing error! Will not affect the video. Use rethrow(draw_error) to see the details!')
end

if(WRITE_VIDEO)
    close(v);
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