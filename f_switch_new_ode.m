IS_SWITCH = 0;

IS_DRAW = 0;
IS_ANIMATE = 0;

global sigma;
sigma=1;
a1 = 3; b1 = 1;   % ellipse constants
a2 = 2; b2 = 3; c2 = 2; % obstacle parameters
xoff = 0; yoff = -1;    % obstacle offset

k1 = 1; k2 = 0.4;    % path-following and reactive vector field gains
k3 = 0.4;            % perturbed reactive vector field gain
l1 = 0.1; l2 = 0.1;             % bump function constants
c = -1.5;              % repulsive bounary \varphi^{-1}(c)
delta = 0.5;            % perturbed reactive boundary 
param={a1, b1, a2, b2, c2, xoff, yoff, k1, k2, k3, l1, l2, c, delta, IS_SWITCH};

%%% singular points in the mixed area; calculated by f_switch_calsingular.nb
singu1 = [-1.08978, -1.89536];
singu2 = [-1.01443, -1.60802];
singu3 = [0.406209, -0.0477639];

%% ode45
opts = odeset('AbsTol',1e-15,'MaxStep',0.1);
tspan = [0 50];
%%% test 1
initial = [2; 0];
[tout1, p1] = ode45(@(t, p) sim(t,p, param), tspan, initial, opts);
[~, e1] = cellfun(@(t, p) sim(t,p, param), num2cell(tout1), num2cell(p1,2), 'UniformOutput',false);
e1 = cell2mat(e1);
%%% test 2
initial = [-1.4; 1];
[tout2, p2] = ode45(@(t, p) sim(t,p, param), tspan, initial, opts);
[~, e2] = cellfun(@(t, p) sim(t,p, param), num2cell(tout2), num2cell(p2,2), 'UniformOutput',false);
e2 = cell2mat(e2);
%%% test 3
initial = [-2.8; -4];
[tout3, p3] = ode45(@(t, p) sim(t,p, param), tspan, initial, opts);
[~, e3] = cellfun(@(t, p) sim(t,p, param), num2cell(tout3), num2cell(p3,2), 'UniformOutput',false);
e3 = cell2mat(e3);

%% plot
tout = tout1(1:800,:);
p = p1(1:800,:);
if(IS_DRAW) 
    f_switch_plot; 
end

%% animate
if(IS_ANIMATE)
    f_switch_animate;
end

%% system equations
function [dpdt, e] = sim(t, p, param)
    global sigma;
    x = p(1); y = p(2);
    [a1, b1, a2, b2, c2, xoff, yoff, k1, k2, k3, l1, l2, c, delta, IS_SWITCH] = param{:};
    E = [0, -1; 1, 0];
    
    e1 = x^2/a1^2 + y^2/b1^2 - 1;   % phi
    n1 = [(2*x)/a1^2; (2*y)/b1^2];
    v1 = E*n1 - k1*e1*n1;

    e2 = a2*(x - xoff)^4 + c2*(y - yoff)^4 - b2*(x - xoff)^2*(y - yoff)^2 - 2; % varphi
    n2 = [4*a2*(x - xoff)^3 - b2*(2*x - 2*xoff)*(y - yoff)^2;  
          4*c2*(y - yoff)^3 - b2*(2*y - 2*yoff)*(x - xoff)^2];
    v2 = E*n2 - k2*e2*n2;
    
    %%% bump fuunctions with different radius mixed vf
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
    vsigma1 = kv1 * v1 / norm(v1) + kv2 * v2 / norm(v2);
    
    if(IS_SWITCH)  
        e3 = e2 - delta;  % varphi'
        n3 = n2;
        vsigma2 = E*n3 - k3*e3*n3;

        epsilon = 0.1;
        if( sigma == 1 && abs(e2 - l2*c/(l1+l2)) <= epsilon )
            sigma = 2;
        end
        %%% use the new switching in v7 paper
        outp_normal = [1, 1];
        if( sigma == 2 && e2 >0 && abs(e1)<0.5 && dot(vsigma1,outp_normal)>0)
            sigma = 1;
        end  

        if(sigma==2)
            dpdt = vsigma2;
        else
            dpdt = vsigma1;
        end
    else
        dpdt = vsigma1;
    end
    
    e = [e1 e2];
end