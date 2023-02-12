WRITE_VIDEO = 0;
WINDOWS = 0;
fontsize = 20;

%% plot the desired paths
figure('Position', [100 100 1024 768]); axis equal; set(gcf,'color','w'); hold on; axis equal; % axis off
syms x y
phi = x^2 / a1^2 + y^2 / b1^2 - 1;   % desired path
fimplicit(phi, 'Color', 'red', 'LineWidth', 2)
% plot(0, 0, 'Marker','.','MarkerSize',10);

plot_region = [-3, 3, -3, 3];
fr = a2 * (x-xoff)^4 - b2*(x-xoff)^2*(y-yoff)^2 + c2*(y-yoff)^4-2;  % reactive boundary 
fimplicit(fr, plot_region, 'Color', 'green', 'LineWidth', 2);

f = fr - delta;  % perturbed reactive boundary 
fimplicit(f, plot_region, 'LineStyle', '--', 'Color', 'green', 'LineWidth', 2);

f = fr - c;  % repulsive booundary
fimplicit(f, plot_region, 'LineStyle', '-.', 'Color', 'black', 'LineWidth', 2);
    
f = fr - l2*c/(l1+l2);   % E
fimplicit(f, plot_region, 'LineStyle', '-.', 'Color', 'blue', 'LineWidth', 2);


%%% plot the singular point!!!
plot(singu1(1), singu1(2), 'Marker','+','MarkerSize', 7, 'Color', 'black', 'LineWidth', 2);
plot(singu2(1), singu2(2), 'Marker','o','MarkerSize', 7, 'Color', 'red', 'LineWidth', 2);
plot(singu3(1), singu3(2), 'Marker','+','MarkerSize', 7, 'Color', 'black', 'LineWidth', 2);

%% plot the vector field
%%% calculate the vectors
E = [0, -1; 1 0];
step = 0.2; x = -3.5 : step : 3.5; y = -4.2 : step : 3.2;
[X, Y] = meshgrid(x,y);
VX3 = zeros(size(X));
VY3 = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        x = X(i,j); y = Y(i,j);
        e1 = x^2/a1^2 + y^2/b1^2 - 1;
        gra1 = [(2*x)/a1^2; (2*y)/b1^2];
        e2 = a2*(x - xoff)^4 + c2*(y - yoff)^4 - b2*(x - xoff)^2*(y - yoff)^2 - 2;
        gra2 = [4*a2*(x - xoff)^3 - b2*(2*x - 2*xoff)*(y - yoff)^2;  
          4*c2*(y - yoff)^3 - b2*(2*y - 2*yoff)*(x - xoff)^2];
        
        vec1 = (E - k1 * e1 * eye(2)) * gra1;
        vec2 = (E - k2 * e2 * eye(2)) * gra2;
        
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

plot(p(1,1),p(1,2),'bo','LineWidth', 1); %%% plot starting position
xlabel('X'); ylabel('Y');
set(findall(gcf,'-property','FontSize'),'FontSize', fontsize);  % set font size
hAxes = gca; 


%% animation
h = animatedline(hAxes, 'Color', 'magenta', 'LineWidth', 2);
hh = animatedline(hAxes, 'Color','blue', 'Marker', '.', 'MarkerSize', 20, 'LineWidth', 2);
% axis([-3 3 -3 3]);

numpoints = size(p,1);
drawnum = 3;
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
        % dis = sprintf("Time (s): %0.3f", e.Time(i+(drawnum-1)));
        dis = sprintf("Time (s): %0.3f", tout(i+(drawnum-1)));
        tt = text(hAxes, -4, 3, dis);

        addpoints(h, p(i:i+(drawnum-1),1), p(i:i+(drawnum-1),2));  % trajectory 
        addpoints(hh, p(i+(drawnum-1),1), p(i+(drawnum-1),2)); % trajectory point

        if(WRITE_VIDEO)
            cdata = print('-f1', '-RGBImage','-r0');
            F(drawnum) = im2frame(cdata);
            writeVideo(v, F(drawnum));
        else
            drawnow;
            pause(0.01) ;
        end

        delete(tt)
        clearpoints(hh)
    end
catch draw_error
    disp('Drawing error! Will not affect the video. Use rethrow(draw_error) to see the details!')
end

if(WRITE_VIDEO)
    close(v);
end