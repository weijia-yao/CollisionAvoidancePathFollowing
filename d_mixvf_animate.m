%% plot the desired paths
figure; axis equal; set(gcf,'color','w'); hold on; grid on; axis equal;
xlabel('X'); ylabel('Y');
hAxes = gca;

r = 2;
a = 1; b = 0.5;
beta = pi/4;

syms x y
phi = x^2 + y^2 -r^2;
fimplicit(phi, 'Color', 'red', 'LineWidth', 2)
plot(0, 0, 'Marker','.','MarkerSize',10);

vphi1 = x^2/a^2 + (y+r)^2/b^2 - 1;
fimplicit(vphi1, 'LineWidth', 2)
plot(0, -r, 'Marker','.','MarkerSize',10);

vphi2 = ((x-r)*cos(beta)+y*sin(beta))^2/a^2 + ((x-r)*sin(beta)-y*cos(beta))^2/b^2 - 1;
fimplicit(vphi2, 'LineWidth', 2)
plot(r, 0, 'Marker','.','MarkerSize',10);

vphi3 = (x+2.6)^2/b^2 + y^2/a^2 - 1;
fimplicit(vphi3, 'LineWidth', 2)
plot(-2.6, 0, 'Marker','.','MarkerSize',10);

vphi4 = (x+1.4)^2/b^2 + y^2/a^2 - 1;
fimplicit(vphi4, 'LineWidth', 2)
plot(-1.4, -0, 'Marker','.','MarkerSize',10);

vphi5 = ((x-0.9)^2+(y-r)^2)*((x+0.9)^2+(y-r)^2) - 0.9;
fimplicit(vphi5, 'LineWidth', 2)
plot(0, 0, 'Marker','.','MarkerSize',10);

%%%% draw the bump function boundary
c1 = -0.72;
vphi1_b = vphi1 - c1;
fimplicit(vphi1_b, 'LineStyle', '-.', 'LineWidth', 2);

c2 = -0.72;
vphi2_b = vphi2 - c2;
fimplicit(vphi2_b, 'LineStyle', '-.', 'LineWidth', 2);

c3 = -0.2;
vphi3_b = vphi3 - c3;
fimplicit(vphi3_b, 'LineStyle', '-.', 'LineWidth', 2);

c4 = -0.2;
vphi4_b = vphi4 - c4;
fimplicit(vphi4_b, 'LineStyle', '-.', 'LineWidth', 2);

c5 = -0.2;
vphi5_b = vphi5 - c5;
fimplicit(vphi5_b, 'LineStyle', '-.', 'LineWidth', 2);

%% plot starting position
% hAxes = gca;
plot(p(1,1),p(1,2),'bo','LineWidth', 1);

h = animatedline(hAxes, 'Color', 'magenta', 'LineWidth', 2);
hh = animatedline(hAxes, 'Color', 'blue', 'Marker', '.', 'MarkerSize', 20);
axis([-3.5 3 -3 3]);

%% animation
numpoints = size(p,1);
drawnum = 3;
pause(0);
for i = 1 : drawnum : numpoints-(drawnum-1)
    dis = sprintf("Time (s): %0.3f", e.Time(i+(drawnum-1)));
    tt = text(hAxes, -0.6, 2.8, dis);
    
    % trajectory 
    addpoints(h, p(i:i+(drawnum-1),1), p(i:i+(drawnum-1),2));
    % trajectory point
    addpoints(hh, p(i+(drawnum-1),1), p(i+(drawnum-1),2));
    
    drawnow
    pause(0.01) 
    
    delete(tt)
    clearpoints(hh)
end