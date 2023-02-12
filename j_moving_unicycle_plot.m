fontsize = 20;
if(IS_UNICYCLE)
    index = [1 273 363 500 550 650];
else
    index = [1 180 228 264 314 414];
end
for i=1:length(index)
    %%% plot desired path, reactive, repulsive boundaries
    figure; hold on; set(gcf,'color','w'); axis equal;
    xlim([-8, 1.5]); ylim([-1.5, 1.5]); 
    syms x y
    phi(x,y) = y - sin(x);   % desired path
    fimplicit(phi, 'Color', 'red', 'LineWidth', 2)
    hold on;
    
    t = tout(index(i))
    f = (x+5-obs_v*t)^2/a^2 + y^2/b^2 - 1;  % reactive boundary 
    reacb = fimplicit(f, 'Color', 'green', 'LineWidth', 2);
    f = (x+5-obs_v*t)^2/a^2 + y^2/b^2 - 1 - c;  % repulsive boundary 
    repulb = fimplicit(f, 'LineStyle', '-.', 'Color', 'black', 'LineWidth', 2);
    f=fill(repulb.XData, repulb.YData,'r');               % obstacle
    f.FaceAlpha=0.5;
    
    if(IS_UNICYCLE)
        %%% plot the traj
        plot(p(1:index(i),1),p(1:index(i),2),'m', 'LineWidth',2);
        plot(p(1,1),p(1,2),'bo','MarkerSize', 8); %%% plot starting position
    
        id = index(i);
        L = 0.4;
        theta = p(id,3);
        point=[p(id,1); p(id,2); 0];
        vel = [dpdt(id,1); dpdt(id,2); 0];
        tript = cal_tri_pts(point, vel, L);
        line_pts = [tript(:,3), tript(:,1), tript(:,2), tript(:,3), tript(:,4), tript(:,1)];
        line(line_pts(1,:), line_pts(2,:), line_pts(3,:), 'LineWidth', 2);
    else
        %%% plot the traj
        plot(p(1:index(i),1),p(1:index(i),2),'black', 'LineWidth',2);
        plot(p(1,1),p(1,2),'bo','MarkerSize', 8); %%% plot starting position
    end
    
    %%% plot the vector field
    E = [0, -1; 1 0];
    step = 0.3; x = -8 : step : 1.5; y = -1.5 : step : 1.5;
    [X, Y] = meshgrid(x,y);
    VX3 = zeros(size(X));
    VY3 = zeros(size(X));
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            e1 = Y(i,j);
            gra1 = [0; 1];
            vec1 = (E - k1 * e1 * eye(2)) * gra1;
            
            e2 = (X(i,j)+5-obs_v*t)^2/a^2 + Y(i,j)^2/b^2 - 1;
            gra2 = [2*(X(i,j)+5-obs_v*t)/a^2; 2*Y(i,j)/b^2];
            pvphipt = -2*obs_v*(X(i,j)+5-obs_v*t)/a^2;
            last_term = 1/norm(gra2)^2 * (-pvphipt - l*e2)*gra2;
            vec2 = (E - k2 * e2 * eye(2)) * gra2 + last_term;
            
            if( e2 < c)
                kv1 = 0;
            else
                if( e2 < 0 )
                    f1 = exp( l1/(c - e2) );
                    f2 = exp( l2/(e2 - 0) );
                    kv1 = f1/(f1+f2);
                else
                    kv1 = 1;
                end
            end
            kv2 = 1 - kv1;
            
            vec = kv1 * vec1 / norm(vec1) + kv2 * vec2 / norm(vec2);
    
            VX3(i,j) = vec(1);
            VY3(i,j) = vec(2);
        end
    end
    
    %%% plot the vector field
    quiver(X,Y,VX3, VY3, 'Color', 'b', 'LineWidth', 0.3);

    xlabel('x'); ylabel('y'); 
end
% set(findall(gcf,'-property','FontSize'),'FontSize', fontsize);  % set font size
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