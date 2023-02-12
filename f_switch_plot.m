fontsize = 20;
%% plot the desired path
%%%% plot two desired paths that are used to generate vfs
figure; hold on; set(gcf,'color','w'); axis equal;
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
replbdy=fimplicit(f, plot_region, 'LineStyle', '-.', 'Color', 'black', 'LineWidth', 2);
f=fill(replbdy.XData, replbdy.YData,'r');               % obstacle
f.FaceAlpha=0.5;
    
f = fr - l2*c/(l1+l2);   % E
fimplicit(f, plot_region, 'LineStyle', '-.', 'Color', 'blue', 'LineWidth', 2);


%%% plot the singular point!!!
plot(singu1(1), singu1(2), 'Marker','+','MarkerSize', 7, 'Color', 'black', 'LineWidth', 2);
plot(singu2(1), singu2(2), 'Marker','o','MarkerSize', 7, 'Color', 'red', 'LineWidth', 2);
plot(singu3(1), singu3(2), 'Marker','+','MarkerSize', 7, 'Color', 'black', 'LineWidth', 2);
    
%% plot the traj
% t=1:length(e.Time);
t=1:length(tout);
plot(p(t,1),p(t,2),'m', 'LineWidth',2);
plot(p(1,1),p(1,2),'bs','LineWidth',2); %%% plot starting position
xlabel('x'); ylabel('y');

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
% xlim([-3.5, 3.5]); ylim([-3.2, 2.1]); 
hold off;

%%
axis([-3.1, 3.1, -2.4, 1.2])
text(-2.5, 0, '$\mathcal{P}$', 'Interpreter', 'latex')
text(0, -0.6, '$\mathcal{Q}$', 'Interpreter', 'latex')
text(-1, 0.3, '$\mathcal{R}$', 'Interpreter', 'latex')
text(-1.5, 0.3, '$\mathcal{R}_{\delta}$', 'Interpreter', 'latex')
text(1.5, -0.4, '$\mathcal{E}$', 'Interpreter', 'latex')
set(findall(gcf,'-property','FontSize'),'FontSize', fontsize);  % set font size
set(gca,'FontSize',10)
