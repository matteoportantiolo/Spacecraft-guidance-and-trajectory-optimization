% Spacecraft Guidance and Navigation
% Assignment # 1
% Author: Matteo Portantiolo

%% SYMS
% clearvars; close all; clc
%
% syms x y vx vy mu rho omega_s t m_s
% 
% r1 = sqrt((x+mu).^2+y.^2);
% r2 = sqrt((x+mu-1).^2+y.^2);
% r3 = sqrt( (x-rho*cos(omega_s*t)).^2 + (y-rho*sin(omega_s*t)).^2 );
% v_sq = vx.^2+vy.^2;
% 
% OM = 1/2*(x.^2+y.^2)+(1-mu)/r1+mu/r2+1/2*mu*(1-mu); %OM=U
% J = 2*OM-v_sq;
% 
% OM4 = OM + m_s/r3 - m_s/(rho)^2*(x*cos(omega_s*t)+y*sin(omega_s*t));
% 
% dOM4dx = x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - (m_s*cos(omega_s*t))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (m_s*(2*x - 2*rho*cos(omega_s*t)))/(2*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2)); % can be substituted with jacobian(OM4,x)
% dOM4dy = y - (m_s*sin(omega_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (m_s*(2*y - 2*rho*sin(omega_s*t)))/(2*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2); % can be substituted with jacobian(OM4,y)
% 
% f = [vx;vy;2*vy+dOM4dx;-2*vx+dOM4dy];
% JJ = jacobian(f,[x y vx vy]); % Assemble the matrix A(t)=dfdx

%% START
clearvars; close all; clc; format long g; cspice_kclear();

% Setting plot options
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter','latex')
set(0, 'defaultAxesFontSize', 28 ,'defaultAxesFontSizeMode', 'manual');

%--------------------------------------------------------------------------
%
% NOTE: Comments in section 2.3: Multiple shooting with N=4 
%
% Select right solver options to obtain the results presented in the report
%  1) Optimized options --> Optimal solution in section 2.3 of the report
%  2) General options --> Intermediate solution in Appendix F
%--------------------------------------------------------------------------

%% EX 2.1: Initial guess

% Parameters
alpha = 0.2*pi;
beta = 1.41;
delta = 4;
t_i = 2;
t_f = t_i+delta;
m_s = 3.28900541e5;
rho = 3.88811143e2;
omega_s = -9.25195985e-1;
omega_em = 2.66186135e-1;
l_em = 3.84405e8;
h_i = 167;
h_f = 100;
DU = 3.84405000e5;
TU = 4.34256461;
VU = 1.02454018;

% Spice
cspice_furnsh('ex02.tm');
fprintf('Number of LSK  kernels: %d\n', cspice_ktotal('lsk'));
fprintf('Number of SPK  kernels: %d\n', cspice_ktotal('spk'));
fprintf('Number of PCK  kernels: %d\n', cspice_ktotal('pck'));
fprintf('Number of CK   kernels: %d\n', cspice_ktotal('ck'));
fprintf('Number of TEXT kernels: %d\n', cspice_ktotal('TEXT'));
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

% Gravitational parameter
GM_e = cspice_bodvrd('Earth','GM',1);
GM_m = cspice_bodvrd('Moon','GM',1);
mu = GM_m/(GM_e+GM_m);

% Earth mean radius
radii_Earth = cspice_bodvrd('Earth','RADII',3);
r_e = (radii_Earth(1)+radii_Earth(3)) / 2;

% Moon radius
radii_Moon = cspice_bodvrd('Moon','RADII',3);
r_m = radii_Moon(1);

% Bodies position
x_E = -mu; y_E = 0;  % Earth
x_M = 1-mu; y_M = 0; % Moon

% Moon state
xxm = [1-mu;0;0;0];

% S/C state
r0 = (r_e+h_i)/DU;
rf = (r_m+h_f)/DU;
v0 = beta*sqrt((1-mu)/r0);
x0 = r0*cos(alpha)-mu;
y0 = r0*sin(alpha);
vx0 = -(v0-r0)*sin(alpha);
vy0 = (v0-r0)*cos(alpha);
xx0 = [x0;y0;vx0;vy0];
fprintf('\nInitial state:\n');
fprintf('%d\n',xx0);

% Propagation
[xxf,PHIf,~, xx, tt] = propagatePBR4BP(t_i,xx0,t_f,mu,rho,omega_s,m_s);

% Departure orbit
[x_dep,y_dep] = circular_orbit(r0, x_E, y_E);

% Arrival orbit
[x_arr,y_arr] = circular_orbit(rf, x_M, y_M);

% Rotating frame plot
figure
hold on
plot(xx(:,1),xx(:,2),'LineWidth',5)
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_M,y_M,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(x_E,y_E, 'Earth', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
text(x_M,y_M, 'Moon', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
xlabel('$x_{rot}$ [-]',FontSize=55)
ylabel('$y_{rot}$ [-]',FontSize=55)
grid on
axis equal
title('Rotating frame','Interpreter','Latex',FontSize = 40)

% Zoomed area
x = xx(:,1);
y = xx(:,2);
x_zoom = [-0.05, 0.02]; % X-range of zoomed area
y_zoom = [-0.035, 0.035]; % Y-range of zoomed area
% Add rectangle to highlight zoomed area in main plot
hold on;
% rectangle('Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
%           'EdgeColor', 'red', 'LineWidth', 1.5);
rectangle('Position', [-0.25, -0.25, 0.5, 0.5], ...
          'EdgeColor', 'red', 'LineWidth', 1.5);
% Create inset axes for zoomed view
ax_inset = axes('Position', [0.3, 0.5, 0.25, 0.25]); % Adjust position and size
box on;
hold on;
% Plot zoomed data in inset axes
idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 5);
hold on
plot(xx0(1),xx0(2),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
plot(x_dep, y_dep, 'LineWidth', 5, 'Color', "#A2142F")
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
xlim(x_zoom);
ylim(y_zoom);
grid on;
% Customize inset plot
set(ax_inset, 'FontSize', 30);
legend(["Transfer leg","Initial position","Departure orbit"],'Location','best','Interpreter','Latex',FontSize = 40)
title('Zoomed View','Interpreter','Latex',FontSize = 40)
xlabel('$x_{rot}$ [-]');
ylabel('$y_{rot}$ [-]');

% Conversion to inertial frame
x_in = NaN(size(xx,1),4);
x_m = NaN(size(xx,1),4);
for i = 1:size(xx,1)
    t = tt(i);
    x_in(i,:) = rot2inertial(xx(i,:),t,mu);
    x_m(i,:) = rot2inertial(xxm,t,mu);
end

% Inertial frame plot
figure
hold on
plot(x_in(:,1),x_in(:,2),'LineWidth',5)
plot(x_m(:,1),x_m(:,2),'LineWidth',5,'LineStyle','--','Color',[0.5 0.5 0.5])
plot(0,0,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_m(end,1),x_m(end,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(0,0, 'Earth', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
text(x_m(end,1),x_m(end,2), 'Moon', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
xlabel('$x_{in}$ [-]',FontSize=55)
ylabel('$y_{in}$ [-]',FontSize=55)
grid on
axis equal
legend(["Transfer leg","Moon motion"],'Location','best','Interpreter','Latex',FontSize = 40)
title('Inertial frame','Interpreter','Latex',FontSize = 40)

% % Zoomed area
% x = x_in(:,1);
% y = x_in(:,2);
% x_zoom = [-0.05, 0.02]; % X-range of zoomed area
% y_zoom = [-0.035, 0.035]; % Y-range of zoomed area
% % Add rectangle to highlight zoomed area in main plot
% hold on;
% % rectangle('Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
% %           'EdgeColor', 'red', 'LineWidth', 1.5);
% rectangle('Position', [-0.2, -0.2, 0.4, 0.4], ...
%           'EdgeColor', 'red', 'LineWidth', 1.5);
% % Create inset axes for zoomed view
% ax_inset = axes('Position', [0.5, 0.6, 0.25, 0.25]); % Adjust position and size
% box on;
% hold on;
% % Plot zoomed data in inset axes
% idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
% plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 2);
% hold on
% plot(x_in0(1),x_in0(2),'o','MarkerEdgeColor','c','MarkerSize',20)
% plot(0,0,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
% xlim(x_zoom);
% ylim(y_zoom);
% grid on;
% % Customize inset plot
% set(ax_inset, 'FontSize', 20);
% 
% title('Zoomed View','Interpreter','Latex',FontSize = 33)
% xlabel('$x_{in}$ [-]');
% ylabel('$y_{in}$ [-]');


%% EX 2.2a: Simple shooting without derivatives

% NO derivatives --> NO STM in integrations and propagations

% Initial guess for state and time
x0tt = [xx0;t_i;t_f];

% Parabolic Earth escape velocity 
v_bound = sqrt(2*(1-mu)/r0);

% Boundaries: alpha=[0,2pi] beta=[1,sqrt(2)] theta=[0,2pi] delta=[0,T] 
% where T=23TU=100days --> adimensional: delta=[0,23TU/TU]=[0,23]
lb = [-mu-r0; -r0; -v_bound; -v_bound; 0; 0];
ub = [-mu+r0; r0; v_bound; v_bound; -2*pi/omega_s; -2*pi/omega_s+23]; % omega_s is negative (relative rotation)

% Linear inequality constraint: (t_i-t_f) < 0
A = [0 0 0 0 1 -1];
b = 0;

% Optimization process
options = optimset('Display','iter','LargeScale','off','Algorithm','active-set','TolCon',1e-10); 
[x0tt_opt, delta_v_opt] = fmincon(@obj_function_noSTM,x0tt,A,b,[],[],lb,ub,@constraint_noSTM,options);

% Results
fprintf('\nOptimal delta velocity without derivatives [-]: %d',delta_v_opt);
fprintf('\nOptimal delta velocity without derivatives [km/s]: %d',delta_v_opt*VU);
fprintf('\nOptimal initial state:\n');
fprintf('%d\n',x0tt_opt);

% Constraints evaluation
[~, c_eq_opt] = constraint_noSTM(x0tt_opt);
fprintf('\nWorst equality constraint without derivatives: %d',max(c_eq_opt));

% Propagation without STM
[~,~, xx, tt] = propagatePBR4BP_noSTM(x0tt_opt(5),x0tt_opt(1:4),x0tt_opt(6),mu,rho,omega_s,m_s);

% Rotating frame plot
figure
hold on
plot(xx(:,1),xx(:,2),'LineWidth',5)
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_M,y_M,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(x_E,y_E, 'Earth', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
text(x_M,y_M, 'Moon', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
xlabel('$x_{rot}$ [-]',FontSize=55)
ylabel('$y_{rot}$ [-]',FontSize=55)
grid on
axis equal
title('Rotating frame','Interpreter','Latex',FontSize = 40)

% Departure zoomed area
x = xx(2:end,1);
y = xx(2:end,2);
x_zoom = [-0.05, +0.02]; % X-range of zoomed area
y_zoom = [-0.035, 0.035]; % Y-range of zoomed area
% Add rectangle to highlight zoomed area in main plot
hold on;
% Create inset axes for zoomed view
ax_inset = axes('Position', [0.28, 0.5, 0.2, 0.2]); % Adjust position and size
box on;
hold on;
% Plot zoomed data in inset axes
idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 5);
hold on
plot(xx(1,1),xx(1,2),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
plot(x_dep, y_dep, 'LineWidth', 5, 'Color', "#A2142F")
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
xlim(x_zoom);
ylim(y_zoom);
grid on;
% Customize inset plot
set(ax_inset, 'FontSize', 20);
legend(["","Initial position","Departure orbit"],'Location','best','Interpreter','Latex',FontSize = 33)
title('Departure Zoom','Interpreter','Latex',FontSize = 40)
xlabel('$x_{rot}$ [-]');
ylabel('$y_{rot}$ [-]');

% Arrival zoomed area
x = xx(:,1);
y = xx(:,2);
x_zoom = [1-mu-0.01, 1-mu+0.02]; % X-range of zoomed area
y_zoom = [-0.015, 0.015]; % Y-range of zoomed area
% Add rectangle to highlight zoomed area in main plot
hold on;
% Create inset axes for zoomed view
ax_inset = axes('Position', [0.55, 0.5, 0.2, 0.2]); % Adjust position and size
box on;
hold on;
% Plot zoomed data in inset axes
idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 5);
hold on
plot(xx(end,1),xx(end,2),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
plot(x_arr, y_arr, 'LineWidth', 5, 'Color', "#EDB120")
plot(x_M,y_M,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
xlim(x_zoom);
ylim(y_zoom);
grid on;
% Customize inset plot
set(ax_inset, 'FontSize', 20);
legend(["","Final position","Arrival orbit"],'Location','best','Interpreter','Latex',FontSize = 33)
title('Arrival Zoom','Interpreter','Latex',FontSize = 40)
xlabel('$x_{rot}$ [-]');
ylabel('$y_{rot}$ [-]');

% Conversion to inertial frame
x_in = NaN(size(xx,1),4);
x_m = NaN(size(xx,1),4);
for i = 1:size(xx,1)
    t = tt(i);
    x_in(i,:) = rot2inertial(xx(i,:),t,mu);
    x_m(i,:) = rot2inertial(xxm,t,mu);
end

% Inertial frame plot
figure
hold on
plot(x_in(:,1),x_in(:,2),'LineWidth',5)
plot(x_m(:,1),x_m(:,2),'r','LineWidth',5,'LineStyle','--','Color',[0.5 0.5 0.5])
plot(0,0,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_m(end,1),x_m(end,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(0,0, 'Earth', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
text(x_m(end,1),x_m(end,2), 'Moon', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
xlabel('$x_{in}$ [-]',FontSize=55)
ylabel('$y_{in}$ [-]',FontSize=55)
grid on
axis equal
legend(["Transfer leg","Moon motion"],'Location','best','Interpreter','Latex',FontSize = 40)
title('Inertial frame','Interpreter','Latex',FontSize = 40)

% % Zoomed area
% x = x_in(:,1);
% y = x_in(:,2);
% % Define zoomed area
% x_zoom = [0.9-0.03, 0.9+0.02]; % X-range of zoomed area
% y_zoom = [0.44-0.025, 0.44+0.025]; % Y-range of zoomed area
% % Add rectangle to highlight zoomed area in main plot
% hold on;
% % rectangle('Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
% %           'EdgeColor', 'red', 'LineWidth', 1.5);
% rectangle('Position', [0-0.1, -0.2-0.1, 0.2, 0.2], ...
%           'EdgeColor', 'red', 'LineWidth', 1.5);
% % Create inset axes for zoomed view
% ax_inset = axes('Position', [0.33, 0.3, 0.25, 0.25]); % Adjust position and size
% box on;
% hold on;
% % Plot zoomed data in inset axes
% idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
% plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 2);
% hold on
% plot(x_in(end,1),x_in(end,2),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
% plot(x_arr, y_arr, 'LineWidth', 2, 'Color', "#A2142F")
% plot(x_m(:,1),x_m(:,2),'r','LineWidth',2,'LineStyle','--','Color',[0.5 0.5 0.5])
% plot(x_m(end,1),x_m(end,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
% xlim(x_zoom);
% ylim(y_zoom);
% grid on;
% % Customize inset plot
% set(ax_inset, 'FontSize', 20);
% title('Zoomed View');
% xlabel('$x_{in}$ [-]');
% ylabel('$y_{in}$ [-]');


%% EX 2.2b: Simple shooting with derivatives

% With derivatives --> STM included in integrations and propagations

% Initial guess for state and time
x0tt = [xx0;t_i;t_f];

% Parabolic Earth escape velocity
v_bound = sqrt(2*(1-mu)/r0);  

% Boundaries: alpha=[0,2pi] beta=[1,sqrt(2)] theta=[0,2pi] delta=[0,T] 
% where T=23TU=100days --> adimensional: delta=[0,23TU/TU]=[0,23]
lb = [-mu-r0; -r0; -v_bound; -v_bound; 0; 0];
ub = [-mu+r0; r0; v_bound; v_bound; -2*pi/omega_s; -2*pi/omega_s+23]; % omega_s is negative (relative rotation)

% Linear inequality constraint: (t_i-t_f) < 0
A = [0 0 0 0 1 -1];
b = 0;

% Optimization process
options = optimset('Display','iter','LargeScale','off','Algorithm','active-set', ...
    'GradObj','on','GradConstr','on','TolCon',1e-10);
[x0tt_opt, delta_v_opt, exitflag,~,~,~,hessian] = fmincon(@obj_function,x0tt,A,b,[],[],lb,ub,@constraint,options);

% Results
fprintf('\nOptimal delta velocity with derivatives [-]: %d',delta_v_opt);
fprintf('\nOptimal delta velocity with derivatives [km/s]: %d',delta_v_opt*VU);
fprintf('\nExit flag: %d\n',exitflag);
fprintf('\nMin eig hessian: %d\n',min(eig(hessian)));
fprintf('\nOptimal initial state:\n');
fprintf('%d\n',x0tt_opt);

% Constraints evaluation
[~, c_eq_opt] = constraint(x0tt_opt);
fprintf('\nWorst equality constraint with derivatives[km/s]: %d\n',max(c_eq_opt));

% propagation with STM
[~,~,~, xx, tt] = propagatePBR4BP(x0tt_opt(5),x0tt_opt(1:4),x0tt_opt(6),mu,rho,omega_s,m_s);

% Rotating frame plot 
figure
hold on
plot(xx(:,1),xx(:,2),'LineWidth',5)
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_M,y_M,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(x_E,y_E, 'Earth', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
text(x_M,y_M, 'Moon', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
xlabel('$x_{rot}$ [-]',FontSize=55)
ylabel('$y_{rot}$ [-]',FontSize=55)
grid on
axis equal
title('Rotating frame','Interpreter','Latex',FontSize = 40)

% Departure zoomed area
x = xx(:,1);
y = xx(:,2);
x_zoom = [-0.05, +0.02]; % X-range of zoomed area
y_zoom = [-0.035, 0.035]; % Y-range of zoomed area
% Add rectangle to highlight zoomed area in main plot
hold on;
% Create inset axes for zoomed view
ax_inset = axes('Position', [0.28, 0.5, 0.2, 0.2]); % Adjust position and size
box on;
hold on;
% Plot zoomed data in inset axes
idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 5);
hold on
plot(xx(1,1),xx(1,2),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
plot(x_dep, y_dep, 'LineWidth', 5, 'Color', "#A2142F")
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
xlim(x_zoom);
ylim(y_zoom);
grid on;
% Customize inset plot
set(ax_inset, 'FontSize', 20);
legend(["","Initial position","Departure orbit"],'Location','best','Interpreter','Latex',FontSize = 33)
title('Departure Zoom','Interpreter','Latex',FontSize = 40)
xlabel('$x_{rot}$ [-]');
ylabel('$y_{rot}$ [-]');

% Arrival zoomed area
x = xx(:,1);
y = xx(:,2);
x_zoom = [1-mu-0.0125, 1-mu+0.0175]; % X-range of zoomed area
y_zoom = [-0.015, 0.015]; % Y-range of zoomed area
% Add rectangle to highlight zoomed area in main plot
hold on;
% Create inset axes for zoomed view
ax_inset = axes('Position', [0.55, 0.5, 0.2, 0.2]); % Adjust position and size
box on;
hold on;
% Plot zoomed data in inset axes
idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 5);
hold on
plot(xx(end,1),xx(end,2),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
plot(x_arr, y_arr, 'LineWidth', 5, 'Color', "#EDB120")
plot(x_M,y_M,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
xlim(x_zoom);
ylim(y_zoom);
grid on;
% Customize inset plot
set(ax_inset, 'FontSize', 20);
legend(["","Final position","Arrival orbit"],'Location','best','Interpreter','Latex',FontSize = 33)
title('Arrival Zoom','Interpreter','Latex',FontSize = 40)
xlabel('$x_{rot}$ [-]');
ylabel('$y_{rot}$ [-]');

% Conversion to inertial
x_in = NaN(size(xx,1),4);
x_m = NaN(size(xx,1),4);
for i = 1:size(xx,1)
    t = tt(i);
    x_in(i,:) = rot2inertial(xx(i,:),t,mu);
    x_m(i,:) = rot2inertial(xxm,t,mu);
end

% Inertial frame plot
figure
hold on
plot(x_in(:,1),x_in(:,2),'LineWidth',5)
plot(x_m(:,1),x_m(:,2),'r','LineWidth',5,'LineStyle','--','Color',[0.5 0.5 0.5])
plot(0,0,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_m(end,1),x_m(end,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(0,0, 'Earth', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
text(x_m(end,1),x_m(end,2), 'Moon', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
xlabel('$x_{in}$ [-]',FontSize=55)
ylabel('$y_{in}$ [-]',FontSize=55)
grid on
axis equal
legend(["Transfer leg","Moon motion"],'Location','best','Interpreter','Latex',FontSize = 40)
title('Inertial frame','Interpreter','Latex',FontSize = 40)

% % Zoomed area
% x = x_in(:,1);
% y = x_in(:,2);
% % Define zoomed area
% x_zoom = [1-0.03, 1+0.02]; % X-range of zoomed area
% y_zoom = [-0.08-0.025, -0.08+0.025]; % Y-range of zoomed area
% % Add rectangle to highlight zoomed area in main plot
% hold on;
% % rectangle('Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
% %           'EdgeColor', 'red', 'LineWidth', 1.5);
% rectangle('Position', [1-0.1, -0.08-0.1, 0.2, 0.2], ...
%           'EdgeColor', 'red', 'LineWidth', 1.5);
% % Create inset axes for zoomed view
% ax_inset = axes('Position', [0.5, 0.6, 0.25, 0.25]); % Adjust position and size
% box on;
% hold on;
% % Plot zoomed data in inset axes
% idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
% plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 2);
% hold on
% plot(x_in(end,1),x_in(end,2),'o','MarkerEdgeColor','m','MarkerSize',20)
% plot(x_m(end,1),x_m(end,2),'.','MarkerEdgeColor','r','MarkerSize',20) 
% plot(x_m(:,1),x_m(:,2),'r','LineWidth',2,'LineStyle','--')
% xlim(x_zoom);
% ylim(y_zoom);
% grid on;
% % Customize inset plot
% set(ax_inset, 'FontSize', 20);
% title('Zoomed View');
% xlabel('$x_{in}$ [-]');
% ylabel('$y_{in}$ [-]');


%% EX 2.3: Multiple shooting with N=4

% Initial guess propagation
Y0 = NaN(18,1);
Y0(1:4) = xx0;
[xx_prop12,~,~, ~, ~]  = propagatePBR4BP(t_i, xx0, t_i+4/3, mu, m_s, rho, omega_s);
Y0(5:8) = xx_prop12;
[xx_prop23,~,~, ~, ~]  = propagatePBR4BP(t_i+4/3, Y0(5:8), t_i+8/3, mu, m_s, rho, omega_s);
Y0(9:12) = xx_prop23;
[xx_prop34,~,~, ~, ~]  = propagatePBR4BP(t_i+8/3, Y0(9:12), t_f, mu, m_s, rho, omega_s);
Y0(13:16) = xx_prop34;
Y0(end-1:end) = [t_i;t_f];

% Linear inequality constraint: A*Y <= b
A = zeros(1, 18);  % Constraint ensures that t0 - tf = 0
A(1, 17) = 1;  
A(1, 18) = -1;
b  = 0;  % b = 0 for the linear inequality constraint

%--------------------------------------------------------------------------
% Select right solver options to obtain the results presented in the report
%  1) Optimized options --> Optimal solution in section 2.3 of the report
%  2) General options --> Intermediate solution in Appendix F

% % 1) Optimized options
% options = optimoptions('fmincon', ...
%     'Display', 'iter-detailed', ...             
%     'Algorithm', 'active-set', ...              
%     'SpecifyObjectiveGradient', true, ...     
%     'SpecifyConstraintGradient', true, ...     
%     'ConstraintTolerance', 1e-10, ...
%     'MaxIterations', 1e4, ...
%     'MaxFunctionEvaluations', 4e4);

% 2) General options
options = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...              
    'Algorithm', 'active-set', ...              
    'SpecifyObjectiveGradient', true, ...     
    'SpecifyConstraintGradient', true, ...     
    'ConstraintTolerance', 1e-10);
%--------------------------------------------------------------------------

% optimization process
[Y0_opt, delta_v_opt_multiple, exitflag,~,~,~,hessian] = fmincon(@obj_fun_multiple,Y0,A,b,[],[],[],[],@constr_multiple,options); %see eig(hessian)>0

% Results
fprintf('\nOptimal delta veocity with multiple shooting [-]: %d\n',delta_v_opt_multiple);
fprintf('\nOptimal delta velocity with multiple shooting [km/s]: %d',delta_v_opt_multiple*VU);
fprintf('\nExit flag: %d\n',exitflag);
fprintf('\nMin eig hessian: %d\n',min(eig(hessian))); %it must be positive 
fprintf('\nOptimal initial state:\n');
fprintf('%d\n',Y0_opt);

% Constraints evaluation
[c_ineq_opt_multiple, c_eq_opt_multiple] = constr_multiple(Y0_opt);
fprintf('\nWorst equality constraint with multiple shooting: %d\n',max(c_eq_opt_multiple));
fprintf('\nWorst inequality constraint with multiple shooting: %d\n',max(c_ineq_opt_multiple));

% Gradients check
[discrepance_Gineq_f,discrepance_Geq_f] = check_grad(Y0_opt);
fprintf('\nWorst equality gradient check with multiple shooting: %d\n',norm(discrepance_Geq_f,inf));
fprintf('\nWorst inequality gradient check with multiple shooting: %d\n',norm(discrepance_Gineq_f,inf));

% Propagation
[~,~,~, xx, tt] = propagatePBR4BP(Y0_opt(17),Y0_opt(1:4),Y0_opt(18),mu,rho,omega_s,m_s);

% Rotating frame plot
figure
hold on
plot(xx(:,1),xx(:,2),'LineWidth',5,'Color',"#0072BD")
plot(Y0_opt(5), Y0_opt(6),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',"#FFFFFF")
plot(Y0_opt(9), Y0_opt(10),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',"#FFFFFF")
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_M,y_M,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(x_E,y_E, 'Earth', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
text(x_M,y_M, 'Moon', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
xlabel('$x_{rot}$ [-]',FontSize=55)
ylabel('$y_{rot}$ [-]',FontSize=55)
grid on
axis equal
title('Rotating frame','Interpreter','Latex',FontSize = 40)
legend(["Transfer leg","Nodes"],'Location','best','Interpreter','Latex',FontSize = 40)

% Departure zoomed area
x = xx(:,1);
y = xx(:,2);
x_zoom = [-0.05, +0.02]; % X-range of zoomed area
y_zoom = [-0.035, 0.035]; % Y-range of zoomed area
% Add rectangle to highlight zoomed area in main plot
hold on;
% Create inset axes for zoomed view
ax_inset = axes('Position', [0.28, 0.3, 0.2, 0.2]); % Adjust position and size
box on;
hold on;
% Plot zoomed data in inset axes
idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 5);
hold on
plot(xx(1,1),xx(1,2),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
plot(x_dep, y_dep, 'LineWidth', 5, 'Color', "#A2142F")
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
xlim(x_zoom);
ylim(y_zoom);
grid on;
% Customize inset plot
set(ax_inset, 'FontSize', 20);
legend(["","Initial position","Departure orbit"],'Location','best','Interpreter','Latex',FontSize = 33)
title('Departure Zoom','Interpreter','Latex',FontSize = 40)
xlabel('$x_{rot}$ [-]');
ylabel('$y_{rot}$ [-]');

% Arrival zoomed area
x = xx(:,1);
y = xx(:,2);
x_zoom = [1-mu-0.0125, 1-mu+0.0175]; % X-range of zoomed area
y_zoom = [-0.015, 0.015]; % Y-range of zoomed area
% Add rectangle to highlight zoomed area in main plot
hold on;
% Create inset axes for zoomed view
ax_inset = axes('Position', [0.55, 0.3, 0.2, 0.2]); % Adjust position and size
box on;
hold on;
% Plot zoomed data in inset axes
idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 5);
hold on
plot(xx(end,1),xx(end,2),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
plot(x_arr, y_arr, 'LineWidth', 5, 'Color', "#EDB120")
plot(x_M,y_M,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
xlim(x_zoom);
ylim(y_zoom);
grid on;
% Customize inset plot
set(ax_inset, 'FontSize', 20);
legend(["","Final position","Arrival orbit"],'Location','best','Interpreter','Latex',FontSize = 33)
title('Arrival Zoom','Interpreter','Latex',FontSize = 40)
xlabel('$x_{rot}$ [-]');
ylabel('$y_{rot}$ [-]');

% Conversion to inertial
Y_in = NaN(size(xx,1),4);
x_m = NaN(size(xx,1),4);
for i = 1:size(xx,1)
    t = tt(i);
    Y_in(i,:) = rot2inertial(xx(i,:),t,mu);
    x_m(i,:) = rot2inertial(xxm,t,mu);
end

% Nodes
t_vect = linspace(Y0_opt(17),Y0_opt(18),4);
node2 = rot2inertial(Y0_opt(5:8),t_vect(2),mu);
node3 = rot2inertial(Y0_opt(9:12),t_vect(3),mu);

% Inertial frame plot
figure
hold on
plot(Y_in(:,1),Y_in(:,2),'LineWidth',5,'Color',"#0072BD")
plot(node2(1), node2(2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',"#FFFFFF")
plot(node3(1), node3(2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',"#FFFFFF")
plot(x_m(:,1),x_m(:,2),'r','LineWidth',5,'LineStyle','--','Color',[0.5 0.5 0.5])
plot(0,0,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_m(end,1),x_m(end,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(0,0, 'Earth', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40); 
text(x_m(end,1),x_m(end,2), 'Moon', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
xlabel('$x_{in}$ [-]',FontSize=55)
ylabel('$y_{in}$ [-]',FontSize=55)
grid on
axis equal
legend(["Transfer leg","Nodes","","Moon motion"],'Location','best','Interpreter','Latex',FontSize = 40)
title('Inertial frame','Interpreter','Latex',FontSize = 40)

% % Zoomed area
% x = x_in(:,1);
% y = x_in(:,2);
% % Define zoomed area
% x_zoom = [-0.03, 0.02]; % X-range of zoomed area
% y_zoom = [-0.025, 0.025]; % Y-range of zoomed area
% % Add rectangle to highlight zoomed area in main plot
% hold on;
% % rectangle('Position', [x_zoom(1), y_zoom(1), diff(x_zoom), diff(y_zoom)], ...
% %           'EdgeColor', 'red', 'LineWidth', 1.5);
% rectangle('Position', [-0.1, -0.1, 0.2, 0.2], ...
%           'EdgeColor', 'red', 'LineWidth', 1.5);
% % Create inset axes for zoomed view
% ax_inset = axes('Position', [0.5, 0.6, 0.25, 0.25]); % Adjust position and size
% box on;
% hold on;
% % Plot zoomed data in inset axes
% idx_zoom = x >= x_zoom(1) & x <= x_zoom(2); % Indices of zoomed data
% plot(x(idx_zoom), y(idx_zoom), 'LineWidth', 2);
% hold on
% plot(x_in0(1),x_in0(2),'o','MarkerEdgeColor','c','MarkerSize',20)
% plot(0,0,'.','MarkerEdgeColor','g','MarkerSize',20) 
% xlim(x_zoom);
% ylim(y_zoom);
% grid on;
% % Customize inset plot
% set(ax_inset, 'FontSize', 20);
% title('Zoomed View');
% xlabel('$x_{in}$ [-]');
% ylabel('$y_{in}$ [-]');

% Simple shooting
delta_ts = (x0tt_opt(6)-x0tt_opt(5)) * TU;
delta_vs = delta_v_opt* VU * 1e3;

% Multiple shooting
delta_tm = (Y0_opt(18)-Y0_opt(17)) * TU;
delta_vm = delta_v_opt_multiple * VU * 1e3;

% Pareto front plot
Im = flip(imread('grafico_pareto_effic.png'));
figure
hold on
image([0,100], [3800,4150], Im)
plot(delta_ts, delta_vs, 'o', 'MarkerSize', 20, 'MarkerEdgeColor',"#4DBEEE", 'MarkerFaceColor',"#4DBEEE")
plot(delta_tm, delta_vm, 'o', 'MarkerSize', 20, 'MarkerEdgeColor',"#EDB120", 'MarkerFaceColor',"#EDB120")
ax = gca; 
ax.FontSize = 25;
xlim([0,100])
ylim([3800,4150])
yticks([0 3800 3850 3900 3950 4000 4050 4100 4150])
xlabel('$\Delta t$ (days)', 'FontSize', 33)
ylabel('$\Delta v$ (m/s)', 'FontSize', 33)
legend(["Simple shooting","Multiple shooting"],'Location','best','Interpreter','Latex',FontSize = 33)


%% EX 2.4: N-body propagation

% Initial epoch 
theta_i = wrapTo2Pi(omega_s*x0tt_opt(5));

% Define list of celestial bodies
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% Initialize propagation data (same as regular n-body)
bodies = nbody_init(labels);

% Select integration frame string (SPICE naming convention)
frame = 'ECLIPJ2000';
center = 'Earth';

% Define initial epoch string
ref_epoch_str = '2024-Sep-28 00:00:00.0000 TDB';
et0 = cspice_str2et(ref_epoch_str);

% Final epoch
final_epoch_str = '2024-Oct-28 00:00:00.0000 TDB';
etf = cspice_str2et(final_epoch_str);

% Time span
et_vect = linspace(et0,etf,10000);
t_vect = datetime(cspice_timout(et_vect, ...
    'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS', 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');

% Theta: relative positions
theta_vect = NaN(size(et_vect));
for i = 1:length(et_vect)

    % Ephemeris
    rv_moon = cspice_spkezr('Moon',et_vect(i),frame,'NONE','EMB');
    rv_sun = cspice_spkezr('Sun',et_vect(i),frame,'NONE','EMB');
    
    % Moon rotating frame 
    x_axis_rot = rv_moon(1:3)/norm(rv_moon(1:3));
    z_axis_rot = cross(rv_moon(1:3),rv_moon(4:6))/norm(cross(rv_moon(1:3),rv_moon(4:6)));
    y_axis_rot = cross(z_axis_rot,x_axis_rot)/norm(cross(z_axis_rot,x_axis_rot));
    Rot = [x_axis_rot,y_axis_rot,z_axis_rot].';

    % Projection of Sun on Moon rotating frame
    r_sun_rot = Rot*rv_sun(1:3);
    theta = wrapTo2Pi(atan2(r_sun_rot(2),r_sun_rot(1)));
    theta_vect(i) = theta;

end

% Filter theta_vect to avoid plotting discontinuities due to 0-->2*pi
theta_diff = abs(diff(theta_vect));
indices = find(theta_diff>1);
theta_vect(indices) = NaN;

% Plot to find guess
figure
plot(t_vect,rad2deg(theta_vect),'LineWidth',4)
hold on
plot([t_vect(1) t_vect(end)],[rad2deg(theta_i) rad2deg(theta_i)],'LineWidth',4,'LineStyle','--')
xlim([t_vect(1) t_vect(end)])
ylim([0 360])
grid on
ylabel('$\theta$ [deg]','Interpreter','latex',FontSize=40)
legend('Angle $\theta$','Initial angle $\theta_i$','Interpreter','latex',FontSize=40)

% Zero finding problem for initial epoch
fun = @(et) initial_epoch(et,theta_i);
et_zero = 7.8203*1e8;
options = optimoptions('fsolve','Display','iter','OptimalityTolerance',1e-13);
et_i = fsolve(fun,et_zero,options);
fprintf('\nInitial epoch:\n');
cspice_et2utc( et_i, 'C', 10 )

% Final istant
delta_t = x0tt_opt(6) - x0tt_opt(5); 
delta_et = delta_t * TU * 86400; 
et_f = et_i + delta_et;
fprintf('\nFinal epoch:\n');
cspice_et2utc( et_f, 'C', 10 )

% State in inertial frame 
x0_in = rot2inertial(x0tt_opt(1:4),x0tt_opt(5),mu);
x0_in = [x0_in(1:2).*DU; 0; x0_in(3:4).*VU; 0]; %3D dimensional 

% Integration
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt, xx] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,center), [et_i et_f], x0_in, options);

% Simple shooting in [km]
x_inn = [x_in(:,1).*DU,x_in(:,2).*DU,zeros(size(x_in,1),1)];

% 3D plot 
figure
hold on
plot3(x_inn(:,1),x_inn(:,2),x_inn(:,3),'LineWidth',5,'Color',"#A2142F") 
plot3(x_inn(end,1),x_inn(end,2),0,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',"#A2142F")
plot3(xx(:,1),xx(:,2),xx(:,3),'LineWidth',5,'Color',"#0072BD") 
plot3(xx(end,1),xx(end,2),xx(end,3),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',"#0072BD")
plot3(0,0,0,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
text(0,0,0, 'Earth', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
xlabel('$x_{in}$ [km]',FontSize=55)
ylabel('$y_{in}$ [km]',FontSize=55)
zlabel('$z_{in}$ [km]',FontSize=55)
grid on
legend(["Simple shooting transfer leg","Simple shooting final position","N-body transfer leg","N-body final position"],'Location','best','Interpreter','Latex',FontSize = 40)
title('3D plot','Interpreter','Latex',FontSize = 40)

% 2D plot 
figure
hold on
plot(x_inn(:,1),x_inn(:,2),'LineWidth',5,'Color',"#A2142F") 
plot(x_inn(end,1),x_inn(end,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',"#A2142F")
plot(xx(:,1),xx(:,2),'LineWidth',5,'Color',"#0072BD") 
plot(xx(end,1),xx(end,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',"#0072BD")
plot(0,0,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
text(0,0, 'Earth', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
xlabel('$x_{in}$ [km]',FontSize=55)
ylabel('$y_{in}$ [km]',FontSize=55)
grid on
title('2D plot','Interpreter','Latex',FontSize = 40)


%% FUNCTIONS

function dxdt = xyPBR4BP_STM(t, xx, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the state derivative and State Transition Matrix (STM) 
%   propagation for the Planar Bicircular Restricted Four-Body Problem 
%   (PBR4BP). The system models the dynamics of a massless particle in the 
%   rotating synodic frame, under the influence of two primaries and a 
%   third perturbing body moving in a circular orbit.
%
% Inputs:
%   t       - Time [non-dimensional]
%   xx      - Extended state vector (20×1):
%               [x; y; vx; vy; Phi(:)]
%             where Phi is the 4×4 State Transition Matrix stored column-wise
%   mu      - Mass parameter, defined as m2 / (m1 + m2) [-]
%   rho     - Distance of the perturbing body from the barycenter [nd]
%   omega_s - Angular velocity of the perturbing body [rad/nd_time]
%   m_s     - Mass of the perturbing body (normalized to total system mass) [-]
%
% Outputs:
%   dxdt - Time derivative of the extended system (20×1):
%             [xdot; ydot; vxdot; vydot; Phidot(:)]
%--------------------------------------------------------------------------

% Extract position and velocity components
x  = xx(1);
y  = xx(2);
vx = xx(3);
vy = xx(4);

% Reshape STM from vector form to 4×4 matrix
Phi = reshape(xx(5:end), 4, 4);

% Partial derivatives of the effective potential (4BP)
dOM4dx = x ...
    - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) ...
    - (m_s*cos(omega_s*t))/rho^2 ...
    + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) ...
    - (m_s*(2*x - 2*rho*cos(omega_s*t))) ...
      /(2*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2));

dOM4dy = y ...
    - (m_s*sin(omega_s*t))/rho^2 ...
    - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) ...
    - (m_s*(2*y - 2*rho*sin(omega_s*t))) ...
      /(2*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2)) ...
    + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

% Jacobian of the equations of motion (4×4 matrix)
dfdx = [ 0, 0, 1, 0
         0, 0, 0, 1
         (mu - 1)/((mu + x)^2 + y^2)^(3/2) ...
           - mu/((mu + x - 1)^2 + y^2)^(3/2) ...
           - m_s/((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2) ...
           + (3*m_s*(2*x - 2*rho*cos(omega_s*t))^2) ...
             /(4*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(5/2)) ...
           + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2)^(5/2)) ...
           - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2)^(5/2)) + 1, ...
           (3*m_s*(2*x - 2*rho*cos(omega_s*t))*(2*y - 2*rho*sin(omega_s*t))) ...
             /(4*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(5/2)) ...
           + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) ...
           - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), ...
           0, 2
         (3*m_s*(2*x - 2*rho*cos(omega_s*t))*(2*y - 2*rho*sin(omega_s*t))) ...
           /(4*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(5/2)) ...
           + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) ...
           - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), ...
           (mu - 1)/((mu + x)^2 + y^2)^(3/2) ...
           - mu/((mu + x - 1)^2 + y^2)^(3/2) ...
           - m_s/((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2) ...
           + (3*m_s*(2*y - 2*rho*sin(omega_s*t))^2) ...
             /(4*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(5/2)) ...
           - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) ...
           + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1, ...
           -2, 0 ];

% STM derivative
Phidot = dfdx * Phi;

% Assemble right-hand side of extended system
dxdt = zeros(20,1);
dxdt(1:2)   = xx(3:4);       % Position derivatives
dxdt(3)     = dOM4dx + 2*vy; % Acceleration in x
dxdt(4)     = dOM4dy - 2*vx; % Acceleration in y
dxdt(5:end) = Phidot(:);     % Flattened STM derivative

end

function dxdt = xyPBR4BP(t, xx, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the state derivative for the Planar Bicircular Restricted 
%   Four-Body Problem (PBR4BP) without the State Transition Matrix (STM).
%   The system models the motion of a massless particle in the synodic 
%   rotating frame, under the influence of two primaries and one additional 
%   perturbing body on a circular orbit.
%
% Inputs:
%   t       - Time [non-dimensional]
%   xx      - State vector (4×1):
%               [x; y; vx; vy]
%             where (x, y) are positions and (vx, vy) are velocities.
%   mu      - Mass parameter, defined as m2 / (m1 + m2) [-]
%   rho     - Distance of the perturbing body from the barycenter [nd]
%   omega_s - Angular velocity of the perturbing body [rad/nd_time]
%   m_s     - Mass of the perturbing body (normalized to system mass) [-]
%
% Outputs:
%   dxdt - Time derivative of the state (4×1):
%             [xdot; ydot; vxdot; vydot]
%--------------------------------------------------------------------------

% Extract state components
x  = xx(1);
y  = xx(2);
vx = xx(3);
vy = xx(4);

% Partial derivatives of the effective potential (4BP)
dOM4dx = x ...
    - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) ...
    - (m_s*cos(omega_s*t))/rho^2 ...
    + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) ...
    - (m_s*(2*x - 2*rho*cos(omega_s*t))) ...
      /(2*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2));

dOM4dy = y ...
    - (m_s*sin(omega_s*t))/rho^2 ...
    - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) ...
    - (m_s*(2*y - 2*rho*sin(omega_s*t))) ...
      /(2*((x - rho*cos(omega_s*t))^2 + (y - rho*sin(omega_s*t))^2)^(3/2)) ...
    + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

% Assemble right-hand side
dxdt = zeros(4,1);
dxdt(1:2) = xx(3:4);     % Position derivatives
dxdt(3)   = dOM4dx + 2*vy; % Acceleration in x
dxdt(4)   = dOM4dy - 2*vx; % Acceleration in y

end

function [xf, PHIf, tf, xx, tt] = propagatePBR4BP(t0, x0, tf, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates the Planar Bicircular Restricted Four-Body Problem (PBR4BP)
%   dynamics together with the State Transition Matrix (STM). Integration
%   is carried out using MATLAB's ODE113 solver.
%
% Inputs:
%   t0       - Initial time [non-dimensional]
%   x0       - Initial state vector (4×1): [x; y; vx; vy]
%   tf       - Final time [non-dimensional]
%   mu       - Mass parameter, defined as m2 / (m1 + m2) [-]
%   rho      - Distance of the perturbing body from the barycenter [nd]
%   omega_s  - Angular velocity of the perturbing body [rad/nd_time]
%   m_s      - Mass of the perturbing body (normalized to system mass) [-]
%
% Outputs:
%   xf    - Final state vector (4×1) at tf
%   PHIf  - Final State Transition Matrix (4×4) at tf
%   tf    - Actual final integration time [non-dimensional]
%   xx    - Full integrated state history (N×20: state + STM)
%   tt    - Time vector corresponding to xx
%--------------------------------------------------------------------------

% Initialize State Transition Matrix at t0
Phi0 = eye(4);

% Extended initial condition: state + STM
x0Phi0 = [x0; Phi0(:)];

% Integration options
options_STM = odeset('reltol', 2.5e-14, 'abstol', 2.5e-14);

% Perform integration with PBR4BP dynamics + STM
[tt, xx] = ode113(@(t, x) xyPBR4BP_STM(t, x, mu, rho, omega_s, m_s), [t0 tf], x0Phi0, options_STM);

% Extract final state and STM
xf   = xx(end, 1:4)';              % Final state vector
PHIf = reshape(xx(end, 5:end),4,4); % Final STM
tf   = tt(end);                     % Final integration time

end

function [xf, tf, xx, tt] = propagatePBR4BP_noSTM(t0, x0, tf, mu, rho, omega_s, m_s)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates the Planar Bicircular Restricted Four-Body Problem (PBR4BP)
%   dynamics without the State Transition Matrix (STM). Integration is
%   performed using MATLAB's ODE113 solver.
%
% Inputs:
%   t0       - Initial time [non-dimensional]
%   x0       - Initial state vector (4×1): [x; y; vx; vy]
%   tf       - Final time [non-dimensional]
%   mu       - Mass parameter, defined as m2 / (m1 + m2) [-]
%   rho      - Distance of the perturbing body from the barycenter [nd]
%   omega_s  - Angular velocity of the perturbing body [rad/nd_time]
%   m_s      - Mass of the perturbing body (normalized to system mass) [-]
%
% Outputs:
%   xf - Final state vector (4×1) at tf
%   tf - Actual final integration time [non-dimensional]
%   xx - Full integrated state history (N×4)
%   tt - Time vector corresponding to xx
%--------------------------------------------------------------------------

% Integration options
options = odeset('reltol', 2.5e-14, 'abstol', 2.5e-14);

% Perform integration with PBR4BP dynamics (no STM)
[tt, xx] = ode113(@(t, x) xyPBR4BP(t, x, mu, rho, omega_s, m_s), [t0 tf], x0, options);

% Extract final state
xf = xx(end, 1:4)';   % Final state vector
tf = tt(end);         % Final integration time

end

function x_in = rot2inertial(x_rot, t, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Converts a state vector from the rotating synodic frame to the inertial
%   barycentric frame in the Circular Restricted Three-Body Problem (CR3BP).
%
% Inputs:
%   x_rot - State vector in rotating frame (4×1): [x; y; vx; vy]
%           Positions and velocities expressed in synodic rotating frame.
%   t     - Time [non-dimensional]
%   mu    - Mass parameter, defined as m2 / (m1 + m2) [-]
%
% Outputs:
%   x_in  - State vector in inertial frame (4×1): [x; y; vx; vy]
%           Positions and velocities expressed in inertial barycentric frame.
%--------------------------------------------------------------------------

% Extract rotating-frame state
x  = x_rot(1);
y  = x_rot(2);
vx = x_rot(3);
vy = x_rot(4);

% Transformation from rotating to inertial frame
x1  = (x + mu) .* cos(t) - y .* sin(t);
y1  = (x + mu) .* sin(t) + y .* cos(t);
vx1 = (vx - y) .* cos(t) - (vy + x + mu) .* sin(t);
vy1 = (vx - y) .* sin(t) + (vy + x + mu) .* cos(t);

% Construct inertial-frame state vector
x_in = [x1; y1; vx1; vy1];

end

function [x, y] = circular_orbit(r, x_center, y_center)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Generates the coordinates of a circular orbit centered at a specified
%   point, useful for plotting reference orbits.
%
% Inputs:
%   r        - Orbital radius of the circular orbit [scalar]
%   x_center - X-coordinate of the circle center (e.g., attractor) [scalar]
%   y_center - Y-coordinate of the circle center (e.g., attractor) [scalar]
%
% Outputs:
%   x - X-coordinates of the circle (500×1 vector)
%   y - Y-coordinates of the circle (500×1 vector)
%--------------------------------------------------------------------------

% Parametric angle samples
theta = linspace(0, 2*pi, 500).';

% Circle parametric equations
x = x_center + r*cos(theta);
y = y_center + r*sin(theta);

end

function K = const()
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Defines and returns a structure containing the main physical and
%   problem-specific constants for the Earth-Moon system in the 
%   Planar Bicircular Restricted Four-Body Problem (PBR4BP).
%
% Inputs:
%   (none)
%
% Outputs:
%   K - Structure containing constants:
%         .mu      - Mass parameter m2 / (m1 + m2) [-]
%         .rho     - Distance of perturbing body from barycenter [nd]
%         .omega_s - Angular velocity of perturbing body [rad/nd_time]
%         .m_s     - Mass of perturbing body (normalized) [-]
%         .r_e     - Earth radius [km]
%         .r_m     - Moon radius [km]
%         .DU      - Distance unit (Earth–Moon distance) [km]
%         .ri      - Initial radial distance (Earth + h_i)/DU [nd]
%         .rf      - Final radial distance (Moon + h_f)/DU [nd]
%         .v_ci    - Initial circular velocity at ri [nd units]
%         .v_cf    - Final circular velocity at rf [nd units]
%--------------------------------------------------------------------------

% Mass and perturbation parameters
K.mu      = 0.012150584269940;   % Mass parameter
K.rho     = 3.88811143e2;        % Distance of perturbing body [nd]
K.omega_s = -9.25195985e-1;      % Angular velocity of perturbing body
K.m_s     = 3.28900541e5;        % Mass of perturbing body

% Altitudes [km]
h_i = 167;   % Initial altitude above Earth
h_f = 100;   % Final altitude above Moon

% Radii and scaling
K.r_e = 6.378136600000000e+03;   % Earth radius [km]
K.r_m = 1.737400000000000e+03;   % Moon radius [km]
K.DU  = 3.84405000e5;            % Distance unit (Earth–Moon distance) [km]

% Normalized distances and circular velocities
K.ri   = (K.r_e + h_i)/K.DU;             % Normalized initial radius
K.rf   = (K.r_m + h_f)/K.DU;             % Normalized final radius
K.v_ci = sqrt((1 - K.mu)/K.ri);          % Initial circular velocity
K.v_cf = sqrt(K.mu/K.rf);                % Final circular velocity

end

function [f_obj, grad_f] = obj_function(x0tt)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Objective function for trajectory optimization in the Planar Bicircular
%   Restricted Four-Body Problem (PBR4BP). The objective minimizes the total
%   ∆v cost at departure and arrival by penalizing deviations from circular
%   velocity magnitudes at the initial and final orbits.
%
% Inputs:
%   x0tt - Decision vector (6×1):
%             [x0; y0; vx0; vy0; t0; tf]
%          where (x0, y0, vx0, vy0) define the initial state in the rotating
%          frame, and (t0, tf) are the initial and final times.
%
% Outputs:
%   f_obj  - Objective function value (scalar): ∆v cost
%   grad_f - Gradient of objective function with respect to x0tt (6×1),
%            returned only if requested.
%--------------------------------------------------------------------------

% Load constants
K = const();

% Extract decision variables
x  = x0tt(1);
y  = x0tt(2);
vx = x0tt(3);
vy = x0tt(4);
t0 = x0tt(5);
tf = x0tt(6);

% Initial state vector
x0 = [x; y; vx; vy];

% Propagate dynamics with STM
[xf, PHIf, ~, ~, ~] = propagatePBR4BP(t0, x0, tf, K.mu, K.rho, K.omega_s, K.m_s);

% Extract final state components
x_f  = xf(1);
y_f  = xf(2);
vx_f = xf(3);
vy_f = xf(4);

% Initial ∆v mismatch (departure circular velocity)
delta_vi = sqrt((vx - y)^2 + (vy + x + K.mu)^2) - K.v_ci;

% Final ∆v mismatch (arrival circular velocity)
delta_vf = sqrt((vx_f - y_f)^2 + (vy_f + x_f + K.mu - 1)^2) - K.v_cf;

% Objective function: sum of ∆v mismatches
f_obj = delta_vi + delta_vf; % NOTE: abs() could be applied if desired

% Gradient (if requested)
if nargout > 1
    grad_f = f_gradient(x0tt, PHIf, xf);
end

end

function f_obj = obj_function_noSTM(x0tt)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Objective function for trajectory optimization in the Planar Bicircular
%   Restricted Four-Body Problem (PBR4BP) without using the State Transition
%   Matrix (STM). The objective minimizes the total ∆v cost at departure and
%   arrival by penalizing deviations from circular velocity magnitudes at
%   the initial and final orbits.
%
% Inputs:
%   x0tt - Decision vector (6×1):
%             [x0; y0; vx0; vy0; t0; tf]
%          where (x0, y0, vx0, vy0) define the initial state in the rotating
%          frame, and (t0, tf) are the initial and final times.
%
% Outputs:
%   f_obj - Objective function value (scalar): ∆v cost
%--------------------------------------------------------------------------

% Load constants
K = const();

% Extract decision variables
x  = x0tt(1);
y  = x0tt(2);
vx = x0tt(3);
vy = x0tt(4);
t0 = x0tt(5);
tf = x0tt(6);

% Initial state vector
x0 = [x; y; vx; vy];

% Propagate dynamics without STM
[xf, ~, ~, ~] = propagatePBR4BP_noSTM(t0, x0, tf, K.mu, K.rho, K.omega_s, K.m_s);

% Extract final state components
x_f  = xf(1);
y_f  = xf(2);
vx_f = xf(3);
vy_f = xf(4);

% Initial ∆v mismatch (departure circular velocity)
delta_vi = sqrt((vx - y)^2 + (vy + x + K.mu)^2) - K.v_ci;

% Final ∆v mismatch (arrival circular velocity)
delta_vf = sqrt((vx_f - y_f)^2 + (vy_f + x_f + K.mu - 1)^2) - K.v_cf;

% Objective function: sum of ∆v mismatches
f_obj = delta_vi + delta_vf; % NOTE: abs() could be applied if desired

end

function grad_f = f_gradient(x0tt, PHIf, xf)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the gradient of the objective function with respect to the
%   decision vector in the Planar Bicircular Restricted Four-Body Problem
%   (PBR4BP). The gradient accounts for the sensitivities of both the
%   initial ∆v mismatch and the propagated final state through the STM.
%
% Inputs:
%   x0tt - Decision vector (6×1):
%             [x0; y0; vx0; vy0; t0; tf]
%          where (x0, y0, vx0, vy0) define the initial state in the rotating
%          frame, and (t0, tf) are the initial and final times.
%   PHIf - Final State Transition Matrix (4×4) at time tf
%   xf   - Final state vector (4×1): [x_f; y_f; vx_f; vy_f]
%
% Outputs:
%   grad_f - Gradient of the objective function (6×1) with respect to x0tt
%--------------------------------------------------------------------------

% Load constants
K = const();

% Extract decision variables
x  = x0tt(1);
y  = x0tt(2);
vx = x0tt(3);
vy = x0tt(4);
t0 = x0tt(5);
tf = x0tt(6);

% Initial state vector
x0 = [x; y; vx; vy];

% Extract final state components
x_f  = xf(1);
y_f  = xf(2);
vx_f = xf(3);
vy_f = xf(4);

% Gradient of initial ∆v mismatch w.r.t. initial state
dvi_dxi = [ (vy + x + K.mu);
            -(vx - y);
             (vx - y);
             (vy + x + K.mu) ] ...
            ./ sqrt((vx - y)^2 + (vy + x + K.mu)^2);

% Gradient of final ∆v mismatch w.r.t. final state
dvf_dxf = [ (vy_f + x_f + K.mu - 1);
            -(vx_f - y_f);
             (vx_f - y_f);
             (vy_f + x_f + K.mu - 1) ] ...
            ./ sqrt((vx_f - y_f)^2 + (vy_f + x_f + K.mu - 1)^2);

% Sensitivities w.r.t. initial state (chain rule via STM)
grad_f(1:4) = dvi_dxi + PHIf' * dvf_dxf;

% Sensitivities w.r.t. initial and final times
dxdt_0 = xyPBR4BP(t0, x0, K.mu, K.rho, K.omega_s, K.m_s);
dxdt_f = xyPBR4BP(tf, xf, K.mu, K.rho, K.omega_s, K.m_s);

grad_f(5) = - dvf_dxf' * PHIf * dxdt_0; % ∂f/∂t0
grad_f(6) =   dvf_dxf' * dxdt_f;        % ∂f/∂tf

end

function [c_ineq, c_eq, grad_c_ineq, grad_c_eq] = constraint(x0tt)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Defines the nonlinear constraints for trajectory optimization in the
%   Planar Bicircular Restricted Four-Body Problem (PBR4BP). The equality
%   constraints enforce that the initial and final states lie on the desired
%   circular orbits (radius and tangency conditions). No inequality
%   constraints are currently imposed.
%
% Inputs:
%   x0tt - Decision vector (6×1):
%             [x0; y0; vx0; vy0; t0; tf]
%          where (x0, y0, vx0, vy0) define the initial state in the rotating
%          frame, and (t0, tf) are the initial and final times.
%
% Outputs:
%   c_ineq     - Inequality constraint vector ([] since none imposed)
%   c_eq       - Equality constraint vector (4×1):
%                  (1) Initial orbit radius condition
%                  (2) Initial tangency condition
%                  (3) Final orbit radius condition
%                  (4) Final tangency condition
%   grad_c_ineq - Gradient of inequality constraints (empty)
%   grad_c_eq   - Gradient of equality constraints (6×4), returned only if requested
%--------------------------------------------------------------------------

% Load constants
K = const();

% Extract decision variables
x  = x0tt(1);
y  = x0tt(2);
vx = x0tt(3);
vy = x0tt(4);
t0 = x0tt(5);
tf = x0tt(6);

% Initial state vector
x0 = [x; y; vx; vy];

% --- Initial constraints (circular orbit at ri) --------------------------
c_eq(1) = (x + K.mu)^2 + y^2 - K.ri^2;                        % Radius
c_eq(2) = (x + K.mu)*(vx - y) + y*(vy + x + K.mu);            % Tangency

% --- Propagation to final state ------------------------------------------
[xf, PHIf, ~, ~, ~] = propagatePBR4BP(t0, x0, tf, K.mu, K.rho, K.omega_s, K.m_s);

x_f  = xf(1);
y_f  = xf(2);
vx_f = xf(3);
vy_f = xf(4);

% --- Final constraints (circular orbit at rf) ----------------------------
c_eq(3) = (x_f + K.mu - 1)^2 + y_f^2 - K.rf^2;                % Radius
c_eq(4) = (x_f + K.mu - 1)*(vx_f - y_f) + y_f*(vy_f + x_f + K.mu - 1); % Tangency

% --- Inequality constraints (none) ---------------------------------------
c_ineq = [];

% --- Gradients -----------------------------------------------------------
if nargout > 2
    grad_c_eq   = c_gradient(x0tt, PHIf, xf); % Gradient of equalities
    grad_c_ineq = [];                         % No inequality constraints
end

end

function [c_ineq, c_eq] = constraint_noSTM(x0tt)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Defines the nonlinear constraints for trajectory optimization in the
%   Planar Bicircular Restricted Four-Body Problem (PBR4BP) without using
%   the State Transition Matrix (STM). The equality constraints enforce
%   that the initial and final states lie on the desired circular orbits
%   (radius and tangency conditions). No inequality constraints are imposed.
%
% Inputs:
%   x0tt - Decision vector (6×1):
%             [x0; y0; vx0; vy0; t0; tf]
%          where (x0, y0, vx0, vy0) define the initial state in the rotating
%          frame, and (t0, tf) are the initial and final times.
%
% Outputs:
%   c_ineq - Inequality constraint vector ([] since none imposed)
%   c_eq   - Equality constraint vector (4×1):
%              (1) Initial orbit radius condition
%              (2) Initial tangency condition
%              (3) Final orbit radius condition
%              (4) Final tangency condition
%--------------------------------------------------------------------------

% Load constants
K = const();

% Extract decision variables
x  = x0tt(1);
y  = x0tt(2);
vx = x0tt(3);
vy = x0tt(4);
t0 = x0tt(5);
tf = x0tt(6);

% Initial state vector
x0 = [x; y; vx; vy];

% --- Initial constraints (circular orbit at ri) --------------------------
c_eq(1) = (x + K.mu)^2 + y^2 - K.ri^2;                        % Radius
c_eq(2) = (x + K.mu)*(vx - y) + y*(vy + x + K.mu);            % Tangency

% --- Propagation to final state ------------------------------------------
[xf, ~, ~, ~] = propagatePBR4BP_noSTM(t0, x0, tf, K.mu, K.rho, K.omega_s, K.m_s);

x_f  = xf(1);
y_f  = xf(2);
vx_f = xf(3);
vy_f = xf(4);

% --- Final constraints (circular orbit at rf) ----------------------------
c_eq(3) = (x_f + K.mu - 1)^2 + y_f^2 - K.rf^2;                % Radius
c_eq(4) = (x_f + K.mu - 1)*(vx_f - y_f) + y_f*(vy_f + x_f + K.mu - 1); % Tangency

% --- Inequality constraints (none) ---------------------------------------
c_ineq = [];

end

function grad_c = c_gradient(x0tt, PHIf, xf)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the gradient of the equality constraints for the trajectory
%   optimization problem in the Planar Bicircular Restricted Four-Body
%   Problem (PBR4BP). The constraints enforce circular orbit conditions at
%   departure and arrival. Gradients are computed w.r.t. the decision
%   vector, including time sensitivities through the STM.
%
% Inputs:
%   x0tt - Decision vector (6×1):
%             [x0; y0; vx0; vy0; t0; tf]
%   PHIf - Final State Transition Matrix (4×4) at time tf
%   xf   - Final state vector (4×1): [x_f; y_f; vx_f; vy_f]
%
% Outputs:
%   grad_c - Gradient matrix of equality constraints (6×4):
%              Each column corresponds to ∇c_i, i = 1...4
%--------------------------------------------------------------------------

% Initialize gradient matrix
grad_c = NaN(4, 6);

% Load constants
K = const();

% Extract decision variables
x  = x0tt(1);
y  = x0tt(2);
vx = x0tt(3);
vy = x0tt(4);
t0 = x0tt(5);
tf = x0tt(6);

% Initial state vector
x0 = [x; y; vx; vy];

% Extract final state components
x_f  = xf(1);
y_f  = xf(2);
vx_f = xf(3);
vy_f = xf(4);

% --- Gradients of initial constraints w.r.t. initial state ---------------
dc1_dxi = [2*(x + K.mu); 2*y; 0; 0];           % Radius condition at t0
dc2_dxi = [vx; vy; x + K.mu; y];               % Tangency condition at t0

% --- Gradients of final constraints w.r.t. final state -------------------
dc3_dxf = [2*(x_f + K.mu - 1); 2*y_f; 0; 0];   % Radius condition at tf
dc4_dxf = [vx_f; vy_f; x_f + K.mu - 1; y_f];   % Tangency condition at tf

% Dynamics at initial and final times
dxdt_0 = xyPBR4BP(t0, x0, K.mu, K.rho, K.omega_s, K.m_s);
dxdt_f = xyPBR4BP(tf, xf, K.mu, K.rho, K.omega_s, K.m_s);

% --- Assemble gradients --------------------------------------------------
grad_c(1,:) = [dc1_dxi', 0, 0];
grad_c(2,:) = [dc2_dxi', 0, 0];
grad_c(3,:) = [(PHIf' * dc3_dxf)', -(PHIf' * dc3_dxf)' * dxdt_0, dc3_dxf' * dxdt_f];
grad_c(4,:) = [(PHIf' * dc4_dxf)', -(PHIf' * dc4_dxf)' * dxdt_0, dc4_dxf' * dxdt_f];

% Return as (6×4) matrix
grad_c = grad_c';

end

function [f_obj, grad_f] = obj_fun_multiple(Y)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Objective function for multi-leg trajectory optimization in the Planar
%   Bicircular Restricted Four-Body Problem (PBR4BP). The objective is the
%   sum of ∆v mismatches at departure and arrival, evaluated from the first
%   and last state vectors in the decision vector Y.
%
% Inputs:
%   Y - Decision vector containing multiple shooting nodes:
%          Y = [x1; ...; xN]
%       where x1 = [x; y; vx; vy] at initial time,
%             xN = [x; y; vx; vy] at final time.
%       (Here only the first and last 4 states are used.)
%
% Outputs:
%   f_obj  - Objective function value (scalar): total ∆v cost
%   grad_f - Gradient of the objective function w.r.t. Y, returned only
%            if requested.
%--------------------------------------------------------------------------

% Load constants
K = const();

% Extract initial and final states from decision vector
x1  = Y(1:4);       % Initial state
x_N = Y(13:16);     % Final state (assumes N ≥ 4 nodes)

% Initial ∆v mismatch
delta_vi = sqrt((x1(3) - x1(2))^2 + (x1(4) + x1(1) + K.mu)^2) - K.v_ci;

% Final ∆v mismatch
delta_vf = sqrt((x_N(3) - x_N(2))^2 + (x_N(4) + x_N(1) + K.mu - 1)^2) - K.v_cf;

% Objective function: sum of ∆v mismatches
f_obj = delta_vi + delta_vf;

% Gradient (if requested)
if nargout > 1
    grad_f = f_grad_multiple(Y);
end

end

function grad_f = f_grad_multiple(Y)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the gradient of the objective function for the multiple-shooting
%   transcription in the Planar Bicircular Restricted Four-Body Problem 
%   (PBR4BP). The gradient accounts only for sensitivities of the departure
%   and arrival ∆v mismatches with respect to the first and last state
%   vectors in Y.
%
% Inputs:
%   Y - Decision vector containing multiple shooting nodes:
%          Y = [x1; y1; vx1; vy1; ...; xN; yN; vxN; vyN]
%       (Here only the first 4 and the last 4 components are used.)
%
% Outputs:
%   grad_f - Gradient of the objective function (18×1 vector).
%            Nonzero entries appear only in:
%               grad_f(1:4)   → derivatives wrt initial state
%               grad_f(13:16) → derivatives wrt final state
%--------------------------------------------------------------------------

% Load constants
K = const();

% Extract initial state components
x1  = Y(1);
y1  = Y(2);
vx1 = Y(3);
vy1 = Y(4);

% Extract final state components
x_N  = Y(13);
y_N  = Y(14);
vx_N = Y(15);
vy_N = Y(16);

% Gradient of initial ∆v mismatch wrt initial state
dv1_dx1 = [ (vy1 + x1 + K.mu);
            -(vx1 - y1);
             (vx1 - y1);
             (vy1 + x1 + K.mu) ] ...
            ./ sqrt((vx1 - y1)^2 + (vy1 + x1 + K.mu)^2);

% Gradient of final ∆v mismatch wrt final state
dvN_dxN = [ (vy_N + x_N + K.mu - 1);
            -(vx_N - y_N);
             (vx_N - y_N);
             (vy_N + x_N + K.mu - 1) ] ...
            ./ sqrt((vx_N - y_N)^2 + (vy_N + x_N + K.mu - 1)^2);

% Initialize gradient vector (assume 18 variables in Y)
grad_f = zeros(18,1);

% Assign partial derivatives to initial and final state entries
grad_f(1:4)   = dv1_dx1;
grad_f(13:16) = dvN_dxN;

end

function [c_ineq, c_eq, grad_c_ineq, grad_c_eq] = constr_multiple(Y)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Defines the nonlinear constraints for a multiple-shooting transcription
%   of the trajectory optimization problem in the Planar Bicircular
%   Restricted Four-Body Problem (PBR4BP). 
%   - Equality constraints enforce trajectory continuity between nodes and 
%     circular orbit conditions at departure and arrival.
%   - Inequality constraints enforce exclusion regions around Earth and 
%     Moon and enforce time ordering.
%
% Inputs:
%   Y - Decision vector containing:
%         [x1; y1; vx1; vy1; ...; xN; yN; vxN; vyN; t1; tN]
%       Specifically:
%         x1..xx_N : state vectors at 4 shooting nodes (4×1 each)
%         t1, tN   : initial and final times
%
% Outputs:
%   c_ineq     - Inequality constraint vector (9×1):
%                  (1–8)  Node exclusion from Earth and Moon radii
%                  (9)    Time ordering condition (t1 ≤ tN)
%   c_eq       - Equality constraint vector (16×1):
%                  (1–12) Trajectory continuity at shooting nodes
%                  (13–14) Initial orbit radius & tangency conditions
%                  (15–16) Final orbit radius & tangency conditions
%   grad_c_ineq - Gradient of inequality constraints (only if requested)
%   grad_c_eq   - Gradient of equality constraints (only if requested)
%--------------------------------------------------------------------------

% Load constants
K = const();

% Extract node states
xx1 = Y(1:4);    x1 = Y(1);   y1 = Y(2);   vx1 = Y(3);   vy1 = Y(4);
xx2 = Y(5:8);    x2 = Y(5);   y2 = Y(6);
xx3 = Y(9:12);   x3 = Y(9);   y3 = Y(10);
xx_N = Y(13:16); x_N = Y(13); y_N = Y(14); vx_N = Y(15); vy_N = Y(16);

% Extract times
t1 = Y(17);
tN = Y(18);

% Time step (equal segments)
h = (tN - t1) / 3;
t2 = t1 + h;
t3 = t2 + h;

% --- Equality constraints: multiple shooting & orbit conditions ----------
c_eq = NaN(16,1);

% Propagate between nodes
[xf2, PHI12, ~, ~, ~] = propagatePBR4BP(t1, xx1, t2, K.mu, K.rho, K.omega_s, K.m_s);
[xf3, PHI23, ~, ~, ~] = propagatePBR4BP(t2, xx2, t3, K.mu, K.rho, K.omega_s, K.m_s);
[xf4, PHI34, ~, ~, ~] = propagatePBR4BP(t3, xx3, tN, K.mu, K.rho, K.omega_s, K.m_s);

xf   = [xf2; xf3; xf4];     % Store propagated states
PHI  = [PHI12; PHI23; PHI34]; % Store STMs (stacked)

% Continuity constraints (matching propagated and decision states)
c_eq(1:4)   = xf2 - xx2;
c_eq(5:8)   = xf3 - xx3;
c_eq(9:12)  = xf4 - xx_N;

% Initial orbit conditions (radius & tangency)
c_eq(13) = (x1 + K.mu)^2 + y1^2 - K.ri^2;
c_eq(14) = (x1 + K.mu)*(vx1 - y1) + y1*(vy1 + x1 + K.mu);

% Final orbit conditions (radius & tangency)
c_eq(15) = (x_N + K.mu - 1)^2 + y_N^2 - K.rf^2;
c_eq(16) = (x_N + K.mu - 1)*(vx_N - y_N) + y_N*(vy_N + x_N + K.mu - 1);

% --- Inequality constraints: planet avoidance & time ordering ------------
c_ineq = NaN(9,1);

% Earth/Moon exclusion at each node (squared distances vs radii)
c_ineq(1) = (K.r_e/K.DU)^2 - (x1 + K.mu)^2     - y1^2;
c_ineq(2) = (K.r_m/K.DU)^2 - (x1 + K.mu - 1)^2 - y1^2;

c_ineq(3) = (K.r_e/K.DU)^2 - (x2 + K.mu)^2     - y2^2;
c_ineq(4) = (K.r_m/K.DU)^2 - (x2 + K.mu - 1)^2 - y2^2;

c_ineq(5) = (K.r_e/K.DU)^2 - (x3 + K.mu)^2     - y3^2;
c_ineq(6) = (K.r_m/K.DU)^2 - (x3 + K.mu - 1)^2 - y3^2;

c_ineq(7) = (K.r_e/K.DU)^2 - (x_N + K.mu)^2     - y_N^2;
c_ineq(8) = (K.r_m/K.DU)^2 - (x_N + K.mu - 1)^2 - y_N^2;

% Time ordering (ensure t1 ≤ tN)
c_ineq(9) = t1 - tN;

% --- Gradients (if requested) --------------------------------------------
if nargout > 2
    grad_c_eq   = c_grad_multiple(Y, xf, PHI);     % Gradient of equalities
    grad_c_ineq = c_grad_ineq_multiple(Y);         % Gradient of inequalities
end

end

function grad_c_eq = c_grad_multiple(Y, xf, PHI)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the Jacobian (gradient matrix) of the equality constraints 
%   used in the multiple-shooting formulation of the Planar Bicircular 
%   Restricted Four-Body Problem (PBR4BP). 
%   The constraints enforce trajectory continuity between nodes as well as 
%   circular orbit conditions at departure and arrival. Sensitivities with 
%   respect to both the state and the time variables are included.
%
% Inputs:
%   Y   - Decision vector (18×1):
%            [x1; y1; vx1; vy1; x2; y2; vx2; vy2; 
%             x3; y3; vx3; vy3; xN; yN; vxN; vyN; t1; tN]
%   xf  - Propagated states at intermediate nodes (12×1):
%            [xf2; xf3; xf4], each block 4×1
%   PHI - Stacked State Transition Matrices (12×4):
%            [PHI12; PHI23; PHI34], each block 4×4
%
% Outputs:
%   grad_c_eq - Gradient matrix of equality constraints (18×16),
%               where each column corresponds to ∇c_i for i = 1...16.
%--------------------------------------------------------------------------

% Initialize gradient matrix
grad_c = zeros(16, 18);

% Load constants
K = const();

% Extract state variables
x1  = Y(1);   y1  = Y(2);   vx1 = Y(3);   vy1 = Y(4);
x_N = Y(13);  y_N = Y(14);  vx_N = Y(15); vy_N = Y(16);

% Extract times
t1 = Y(17);
tN = Y(18);
h  = (tN - t1) / 3;
t2 = t1 + h;
t3 = t2 + h;

% Extract propagated states
xf2 = xf(1:4);
xf3 = xf(5:8);
xf4 = xf(9:12);

% Extract block STMs
PHI12 = PHI(1:4, 1:4);
PHI23 = PHI(5:8, 1:4);
PHI34 = PHI(9:12, 1:4);

% --- Continuity constraints wrt states ----------------------------------
grad_c(1:4,1:4)    = PHI12;
grad_c(5:8,5:8)    = PHI23;
grad_c(9:12,9:12)  = PHI34;
grad_c(1:4,5:8)    = -eye(4);
grad_c(5:8,9:12)   = -eye(4);
grad_c(9:12,13:16) = -eye(4);

% --- Boundary orbit constraints wrt states ------------------------------
dc1_dx1 = [2*(x1 + K.mu), 2*y1, 0, 0];
dc2_dx1 = [vx1, vy1, x1 + K.mu, y1];

dc3_dxN = [2*(x_N + K.mu - 1), 2*y_N, 0, 0];
dc4_dxN = [vx_N, vy_N, x_N + K.mu - 1, y_N];

grad_c(13,1:4)   = dc1_dx1;
grad_c(14,1:4)   = dc2_dx1;
grad_c(15,13:16) = dc3_dxN;
grad_c(16,13:16) = dc4_dxN;

% --- Continuity constraints wrt times -----------------------------------
dxdt_1  = xyPBR4BP(t1, Y(1:4),   K.mu, K.rho, K.omega_s, K.m_s);
dxdt_12 = xyPBR4BP(t2, xf2,      K.mu, K.rho, K.omega_s, K.m_s);

dxdt_2  = xyPBR4BP(t2, Y(5:8),   K.mu, K.rho, K.omega_s, K.m_s);
dxdt_23 = xyPBR4BP(t3, xf3,      K.mu, K.rho, K.omega_s, K.m_s);

dxdt_3  = xyPBR4BP(t3, Y(9:12),  K.mu, K.rho, K.omega_s, K.m_s);
dxdt_34 = xyPBR4BP(tN, xf4,      K.mu, K.rho, K.omega_s, K.m_s);

grad_c(1:4,17)   = -PHI12*dxdt_1  + (2/3)*dxdt_12;
grad_c(5:8,17)   = -(2/3)*PHI23*dxdt_2 + (1/3)*dxdt_23;
grad_c(9:12,17)  = -(1/3)*PHI34*dxdt_3;

grad_c(1:4,18)   = (1/3)*dxdt_12;
grad_c(5:8,18)   = -(1/3)*PHI23*dxdt_2 + (2/3)*dxdt_23;
grad_c(9:12,18)  = -(2/3)*PHI34*dxdt_3 + dxdt_34;

% --- Boundary orbit constraints wrt times -------------------------------
% (all zeros, already satisfied by formulation)

% Return as (18×16) matrix
grad_c_eq = grad_c';

end

function grad_c_ineq = c_grad_ineq_multiple(Y)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the Jacobian (gradient matrix) of the inequality constraints 
%   for the multiple-shooting formulation in the Planar Bicircular 
%   Restricted Four-Body Problem (PBR4BP). 
%   The constraints enforce minimum distance from Earth and Moon at each 
%   node and ensure correct time ordering (t1 ≤ tN).
%
% Inputs:
%   Y - Decision vector (18×1):
%          [x1; y1; vx1; vy1; x2; y2; vx2; vy2;
%           x3; y3; vx3; vy3; xN; yN; vxN; vyN; t1; tN]
%
% Outputs:
%   grad_c_ineq - Gradient matrix of inequality constraints (18×9).
%                 Each column corresponds to ∇c_i for i = 1...9.
%--------------------------------------------------------------------------

% Initialize gradient matrix
grad_c_ineq = zeros(9,18);

% Load constants
K = const();

% Extract positions at nodes
x1 = Y(1);   y1 = Y(2);
x2 = Y(5);   y2 = Y(6);
x3 = Y(9);   y3 = Y(10);
xN = Y(13);  yN = Y(14);

% --- Earth/Moon exclusion constraints at nodes ---------------------------
grad_c_ineq(1:2,1:4)   = [ -2*(x1+K.mu),    -2*y1, 0, 0;
                           -2*(x1+K.mu-1), -2*y1, 0, 0 ];

grad_c_ineq(3:4,5:8)   = [ -2*(x2+K.mu),    -2*y2, 0, 0;
                           -2*(x2+K.mu-1), -2*y2, 0, 0 ];

grad_c_ineq(5:6,9:12)  = [ -2*(x3+K.mu),    -2*y3, 0, 0;
                           -2*(x3+K.mu-1), -2*y3, 0, 0 ];

grad_c_ineq(7:8,13:16) = [ -2*(xN+K.mu),    -2*yN, 0, 0;
                           -2*(xN+K.mu-1), -2*yN, 0, 0 ];

% --- Time ordering constraint (t1 - tN ≤ 0) ------------------------------
grad_c_ineq(9,17) =  1;
grad_c_ineq(9,18) = -1;

% Return as (18×9) matrix
grad_c_ineq = grad_c_ineq';

end

function [discrepance_Gineq, discrepance_Geq] = check_grad(Y0)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Verifies the correctness of analytical constraint gradients by
%   comparing them with numerical gradients obtained via centered finite
%   differences. The function computes discrepancies for both inequality
%   and equality constraints in the multiple-shooting transcription of the
%   PBR4BP optimization problem.
%
% Inputs:
%   Y0 - Decision vector (18×1):
%          [x1; y1; vx1; vy1; ...; xN; yN; vxN; vyN; t1; tN]
%
% Outputs:
%   discrepance_Gineq - Difference between numerical and analytical
%                       inequality constraint gradients (18×9)
%   discrepance_Geq   - Difference between numerical and analytical
%                       equality constraint gradients (18×16)
%--------------------------------------------------------------------------

% Reference decision vector
test = Y0;

% Initialize storage for finite-difference gradients
dc   = NaN(9,18);    % Inequality constraints
dceq = NaN(16,18);   % Equality constraints

% --- Finite-difference gradient estimation -------------------------------
for ii = 1:18
    variation = sqrt(eps) * max(test(ii), 1);
    variation_vect = zeros(length(test),1);
    variation_vect(ii) = variation;

    % Evaluate constraints at perturbed states
    [temp_c_p, temp_ceq_p] = constr_multiple(test + variation_vect);
    [temp_c_m, temp_ceq_m] = constr_multiple(test - variation_vect);

    % Centered finite difference
    dc(:,ii)   = (temp_c_p   - temp_c_m)   ./ (2*variation);
    dceq(:,ii) = (temp_ceq_p - temp_ceq_m) ./ (2*variation);
end

% --- Analytical gradients from constraint functions ----------------------
[~, ~, G, Geq] = constr_multiple(test);

% --- Discrepancies -------------------------------------------------------
discrepance_Gineq = dc'   - G;    % Inequality constraints
discrepance_Geq   = dceq' - Geq;  % Equality constraints

end

function [bodies] = nbody_init(labels)
%NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (PCK kernels)
%

% Initialize output
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)
    % Store body label
    bodies{i}.name = labels{i};
    % Store body gravitational constant
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
end

end

function [dxdt] = nbody_shift_rhs(t, x, bodies, frame, center)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo (based on A. Morselli, Politecnico di Milano)
% Description:
%   Evaluates the right-hand side of a Newtonian N-body propagator.
%   The reference frame must be J2000 or ECLIPJ2000, and the integration
%   center must be one of the specified bodies. Only Newtonian gravitational
%   accelerations are considered.
%
% Inputs:
%   t      - Ephemeris time (seconds past J2000 TDB) [scalar]
%   x      - Cartesian state vector wrt center body (6×1) [km, km/s]
%   bodies - Cell array of structs created with nbody_init
%   frame  - Reference frame ('J2000' or 'ECLIPJ2000')
%   center - Label of central body (string, must be in bodies)
%
% Outputs:
%   dxdt   - Time derivative of state (6×1):
%              [vx; vy; vz; ax; ay; az]
%--------------------------------------------------------------------------

% Validate inputs
if not( strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000') )
    msg = 'Invalid integration reference frame, select either J2000 or ECLIPJ2000';
    error(msg);
end
if not( any(strcmp(center,{cell2mat(bodies).name})) )
    msg = 'Invalid center selected, select one of the bodies';
    error(msg);
end

% Initialize RHS
dxdt = zeros(6,1);

% Kinematics
dxdt(1:3) = x(4:6);

% Position wrt center
rr_center_obj = x(1:3);

% GM of central body
GM_center = cspice_bodvrd(center,'GM',3);

% Central body acceleration
dist2 = dot(rr_center_obj, rr_center_obj);
dist  = sqrt(dist2);
dxdt(4:6) = dxdt(4:6) - GM_center * rr_center_obj /(dist*dist2);

% Contributions from perturbing bodies
for i = 1:length(bodies)
    if strcmp(bodies{i}.name, center)
        continue
    end
    % Position and velocity of perturbing body wrt center
    rv_center_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', center);
    rho = rv_center_body(1:3);
    r   = rr_center_obj;
    d   = r - rho;

    % Non-inertial correction terms
    q = (dot(r, (r - 2*rho))) / (dot(rho,rho));
    f = (q*(3+3*q+q^2)) / (1 + (1+q)^(3/2));

    % Perturbation acceleration
    aa_grav = - bodies{i}.GM .* (r + rho.*f) / norm(d)^3;

    dxdt(4:6) = dxdt(4:6) + aa_grav;
end

end

function [delta_theta] = initial_epoch(et, theta_target)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the angular offset between the Sun and a target direction
%   in the Moon-centered rotating frame at a given epoch. The Moon's
%   rotating frame is defined by its instantaneous position and velocity
%   with respect to the Earth-Moon barycenter (EMB).
%
% Inputs:
%   et           - Ephemeris time (seconds past J2000 TDB)
%   theta_target - Desired target angle [rad]
%
% Outputs:
%   delta_theta  - Angular difference between the Sun's projection in the
%                  Moon-rotating frame and the target angle [rad]
%--------------------------------------------------------------------------

frame  = 'ECLIPJ2000';
center = 'EMB';

% Retrieve Moon and Sun states wrt Earth-Moon barycenter
rv_moon = cspice_spkezr('Moon', et, frame, 'NONE', center);
rv_sun  = cspice_spkezr('Sun',  et, frame, 'NONE', center);

% --- Define Moon rotating frame ------------------------------------------
x_axis_rot = rv_moon(1:3) / norm(rv_moon(1:3));                   % Moon radial direction
z_axis_rot = cross(rv_moon(1:3), rv_moon(4:6));                   % Normal to orbital plane
z_axis_rot = z_axis_rot / norm(z_axis_rot);
y_axis_rot = cross(z_axis_rot, x_axis_rot);                       % In-plane orthogonal
y_axis_rot = y_axis_rot / norm(y_axis_rot);

% Rotation matrix inertial → rotating frame
Rot = [x_axis_rot, y_axis_rot, z_axis_rot].';

% --- Sun direction in Moon rotating frame -------------------------------
r_sun_rot = Rot * rv_sun(1:3);

% Compute Sun angle in rotating frame
theta = wrapTo2Pi(atan2(r_sun_rot(2), r_sun_rot(1)));

% Offset from target angle
delta_theta = theta - theta_target;

end

