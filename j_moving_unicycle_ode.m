% using the 2D Dubin's car model; moving obstalces

IS_DRAW = 0;
IS_ANIMATE = 1;
IS_UNICYCLE = 1;

global last_theta;
a = 2; b = 1;   % ellipse constants
k1 = 1; k2 = 1;  % path-following and reactive vector field gains
l1 = 0.1; l2 = 0.1;             % bump function constants
c = -0.72;          % repulsive bounary \varphi^{-1}(c)
obs_v = 0.5;            % speed of the obstacle
l = 10;              % gain for the time-varying term
speed = 1;          % unicycle speed
k_theta = 3;        % gain for angular speed control
param={a, b, obs_v, k1, k2, l1, l2, c, l, speed, k_theta, IS_UNICYCLE};

%% ode45
opts = odeset('AbsTol',1e-15,'MaxStep',0.1);
tspan = [0 50];
if(IS_UNICYCLE)
    initial = [0.8; -0.6; 0];
    last_theta = initial(end);
else
    initial = [0.8; -0.6];
end
[tout1, p1] = ode45(@(t, p) sim(t,p, param), tspan, initial, opts);
[dpdt, e1] = cellfun(@(t, p) sim(t,p, param), num2cell(tout1), num2cell(p1,2), 'UniformOutput',false);
e1 = cell2mat(e1);
dpdt = cell2mat(dpdt')';


%% plot
tout = tout1;
p = p1;
if(IS_DRAW) 
    j_moving_unicycle_plot; 
end

%% animate
if(IS_ANIMATE)
    j_moving_unicycle_animate;
end

%% system equations
function [dpdt, e] = sim(t, p, param)
    global last_theta;
    [a, b, obs_v, k1, k2, l1, l2, c, l, speed, k_theta, IS_UNICYCLE] = param{:};
    E = [0, -1; 1, 0];
    if(IS_UNICYCLE)
        x = p(1); y = p(2); theta = p(3);
    else
        x = p(1); y = p(2);
    end
    
    e1 = y - sin(x);   % phi
    n1 = [-cos(x); 1];
    vf1 = E*n1 - k1*e1*n1;

    e2 = (x+5-obs_v*t)^2/a^2 + y^2/b^2 - 1; % varphi
    n2 = [2*(x+5-obs_v*t)/a^2; 2*y/b^2];
    pvphipt = -2*obs_v*(x+5-obs_v*t)/a^2;
    vflast_term = 1/norm(n2)^2 * (-pvphipt - l*e2)*n2;
    vf2 = E*n2 - k2*e2*n2 + vflast_term;
    
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
    vfcom = kv1 * vf1 / norm(vf1) + kv2 * vf2 / norm(vf2);

    if(IS_UNICYCLE)
        xidot=[speed*cos(last_theta); 
               speed*sin(last_theta)];  
        if( e2 < c)
            Jv = [ - (k2*(2*x - 2*obs_v*t + 10)^2)/a^4 - (2*k2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2,                   - 2/b^2 - (2*k2*y*(2*x - 2*obs_v*t + 10))/(a^2*b^2);
                                            2/a^2 - (2*k2*y*(2*x - 2*obs_v*t + 10))/(a^2*b^2), - (4*k2*y^2)/b^4 - (2*k2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2];
            vfnorm = vfcom/norm(vfcom);       % actually, here norm(vfcom)=1                
            dotthetad = -transpose(vfnorm)*E*(eye(2)-vfnorm*transpose(vfnorm))*Jv*xidot/norm(vf2);  % note, norm(vf2)
        else
            if( e2 < 0 )
                Jv = [ (((k2*(2*x - 2*obs_v*t + 10)^2)/a^4 + (2*k2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*(exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))) - 1))/(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(1/2) - (((exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*((l1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(2*x - 2*obs_v*t + 10))/(a^2*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2) - (l2*exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))*(2*x - 2*obs_v*t + 10))/(a^2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)^2)))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))^2 - (l1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(2*x - 2*obs_v*t + 10))/(a^2*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2))*((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2))/(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(1/2) - (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(k1*cos(x)^2 + k1*sin(x)*(y - sin(x))))/((abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))) - ((2*abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*sign((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*((k2*(2*x - 2*obs_v*t + 10)^2)/a^4 + (2*k2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2) + 2*abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*sign((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*(2/a^2 - (2*k2*y*(2*x - 2*obs_v*t + 10))/(a^2*b^2)))*(exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))) - 1)*((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2))/(2*(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(3/2)) + (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(k1*cos(x)*(y - sin(x)) - 1)*(2*abs(cos(x) + k1*(y - sin(x)))*sign(cos(x) + k1*(y - sin(x)))*(sin(x) + k1*cos(x)) - 2*abs(1 - k1*cos(x)*(y - sin(x)))*sign(1 - k1*cos(x)*(y - sin(x)))*(k1*cos(x)^2 + k1*sin(x)*(y - sin(x)))))/(2*(abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(3/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))) - (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(k1*cos(x)*(y - sin(x)) - 1)*((l1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(2*x - 2*obs_v*t + 10))/(a^2*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2) - (l2*exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))*(2*x - 2*obs_v*t + 10))/(a^2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)^2)))/((abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))^2) + (l1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(k1*cos(x)*(y - sin(x)) - 1)*(2*x - 2*obs_v*t + 10))/(a^2*(abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2),  ((2/b^2 + (2*k2*y*(2*x - 2*obs_v*t + 10))/(a^2*b^2))*(exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))) - 1))/(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(1/2) - (((exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*((2*l1*y*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))/(b^2*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2) - (2*l2*y*exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)))/(b^2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)^2)))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))^2 - (2*l1*y*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))/(b^2*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2))*((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2))/(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(1/2) + ((2*abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*sign((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*((4*k2*y^2)/b^4 + (2*k2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2) - 2*abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*sign((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*(2/b^2 + (2*k2*y*(2*x - 2*obs_v*t + 10))/(a^2*b^2)))*(exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))) - 1)*((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2))/(2*(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(3/2)) - (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(k1*cos(x)*(y - sin(x)) - 1)*((2*l1*y*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))/(b^2*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2) - (2*l2*y*exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)))/(b^2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)^2)))/((abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))^2) - (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(2*k1*abs(cos(x) + k1*(y - sin(x)))*sign(cos(x) + k1*(y - sin(x))) - 2*k1*abs(1 - k1*cos(x)*(y - sin(x)))*sign(1 - k1*cos(x)*(y - sin(x)))*cos(x))*(k1*cos(x)*(y - sin(x)) - 1))/(2*(abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(3/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))) + (k1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*cos(x))/((abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))) + (2*l1*y*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(k1*cos(x)*(y - sin(x)) - 1))/(b^2*(abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2);
                                                                      (((exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*((l1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(2*x - 2*obs_v*t + 10))/(a^2*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2) - (l2*exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))*(2*x - 2*obs_v*t + 10))/(a^2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)^2)))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))^2 - (l1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(2*x - 2*obs_v*t + 10))/(a^2*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2))*((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2))/(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(1/2) - ((2/a^2 - (2*k2*y*(2*x - 2*obs_v*t + 10))/(a^2*b^2))*(exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))) - 1))/(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(1/2) + ((2*abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*sign((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*((k2*(2*x - 2*obs_v*t + 10)^2)/a^4 + (2*k2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2) + 2*abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*sign((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*(2/a^2 - (2*k2*y*(2*x - 2*obs_v*t + 10))/(a^2*b^2)))*(exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))) - 1)*((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2))/(2*(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(3/2)) + (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(sin(x) + k1*cos(x)))/((abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))) - (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(cos(x) + k1*(y - sin(x)))*(2*abs(cos(x) + k1*(y - sin(x)))*sign(cos(x) + k1*(y - sin(x)))*(sin(x) + k1*cos(x)) - 2*abs(1 - k1*cos(x)*(y - sin(x)))*sign(1 - k1*cos(x)*(y - sin(x)))*(k1*cos(x)^2 + k1*sin(x)*(y - sin(x)))))/(2*(abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(3/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))) + (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(cos(x) + k1*(y - sin(x)))*((l1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(2*x - 2*obs_v*t + 10))/(a^2*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2) - (l2*exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))*(2*x - 2*obs_v*t + 10))/(a^2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)^2)))/((abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))^2) - (l1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(cos(x) + k1*(y - sin(x)))*(2*x - 2*obs_v*t + 10))/(a^2*(abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2), (((4*k2*y^2)/b^4 + (2*k2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*(exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))) - 1))/(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(1/2) + (((exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*((2*l1*y*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))/(b^2*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2) - (2*l2*y*exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)))/(b^2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)^2)))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))^2 - (2*l1*y*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))/(b^2*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2))*((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2))/(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(1/2) - (k1*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))/((abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))) - ((2*abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*sign((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)*((4*k2*y^2)/b^4 + (2*k2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2) - 2*abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*sign((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)*(2/b^2 + (2*k2*y*(2*x - 2*obs_v*t + 10))/(a^2*b^2)))*(exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))/(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))) - 1)*((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2))/(2*(abs((2*x - 2*obs_v*t + 10)/a^2 - (2*k2*y*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/b^2)^2 + abs((2*y)/b^2 + (k2*(2*x - 2*obs_v*t + 10)*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1))/a^2)^2)^(3/2)) + (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(cos(x) + k1*(y - sin(x)))*((2*l1*y*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))/(b^2*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2) - (2*l2*y*exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)))/(b^2*((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)^2)))/((abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))^2) + (exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(2*k1*abs(cos(x) + k1*(y - sin(x)))*sign(cos(x) + k1*(y - sin(x))) - 2*k1*abs(1 - k1*cos(x)*(y - sin(x)))*sign(1 - k1*cos(x)*(y - sin(x)))*cos(x))*(cos(x) + k1*(y - sin(x))))/(2*(abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(3/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))) - (2*l1*y*exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1))*(cos(x) + k1*(y - sin(x))))/(b^2*(abs(cos(x) + k1*(y - sin(x)))^2 + abs(1 - k1*cos(x)*(y - sin(x)))^2)^(1/2)*(exp(l2/((x - obs_v*t + 5)^2/a^2 + y^2/b^2 - 1)) + exp(l1/(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)))*(c - (x - obs_v*t + 5)^2/a^2 - y^2/b^2 + 1)^2)];
                vfnorm = vfcom/norm(vfcom);                     
                norm(vfcom) 
                dotthetad = -transpose(vfnorm)*E*(eye(2)-vfnorm*transpose(vfnorm))*Jv*xidot/norm(vfcom);
            else
                Jv = [  - k1*cos(x)^2 - k1*sin(x)*(y - sin(x)), k1*cos(x);
                         sin(x) + k1*cos(x),       -k1];
                vfnorm = vfcom/norm(vfcom);       % actually, here norm(vfcom)=1                
                dotthetad = -transpose(vfnorm)*E*(eye(2)-vfnorm*transpose(vfnorm))*Jv*xidot/norm(vf1);  % note, norm(vf1)
            end
        end


        hp = [cos(theta); sin(theta)];
        theta_u = dotthetad - k_theta*transpose(hp)*E*vfcom/norm(vfcom);

        dotx = speed*cos(theta);
        doty = speed*sin(theta);
        dottheta = theta_u;
        dpdt = [dotx; doty; dottheta];
        
        last_theta = theta;
    else
        dpdt = vfcom;
    end
    e = [e1 e2];    
end
