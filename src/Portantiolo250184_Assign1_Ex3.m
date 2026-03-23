% Spacecraft Guidance and Navigation
% Assignment # 1
% Author: Matteo Portantiolo

%% START

clearvars; close all; clc; format long g;

% Setting plot options
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter','latex')
set(0, 'defaultAxesFontSize', 28 ,'defaultAxesFontSizeMode', 'manual');

%--------------------------------------------------------------------------
%
% NOTE: Comments in section 3.4: Solve the problem
%
% The while loop to find X0_best (lambda0 + tf) is commented out
% The solution X0_best of the while loop is used afterward 
%--------------------------------------------------------------------------

%% EX 3.1: Debris spatial density

% Data
hi = 800;
hf = 1000; 
delta_i = deg2rad(0.75); 
Re = 6378.1366;
mu = 398600.435; 
rho0 = 750 + Re; 
k1 = 1e-5; 
k2 = 1e-4; 
m0 = 1000; 
Tmax = 3; 
Is = 3120; 
DU = 7178.1366; 
MU = m0; 
g0 = 0.00980665; % [km/s^2] Earth's standard gravity at sea level

% Debry density function
q_fun = @(rho) k1./(k2 + ((rho - rho0)./DU).^2);
h_vect = linspace(hi-100,hf+100,1000);
rho_vect = Re+h_vect;
q_vect = q_fun(rho_vect);

% Debry density plot
figure
plot(h_vect,q_vect,'LineWidth',4)
xlabel('Altitude [km]',FontSize=40)
ylabel('Debris densty [$DU^-3$]',FontSize=40)
grid on

% To see region of debris traffic --> larger scale
figure
h_vect2 = linspace(hi-1000,hf+1000,1000);
rho_vect2 = Re+h_vect2;
plot(h_vect2,q_fun(rho_vect2),'LineWidth',4)
xlabel('Altitude [km]',FontSize=40)
ylabel('Debris densty [$DU^-3$]',FontSize=40)
grid on 

% Initial state in J2000 frame
ri = hi + Re;
vi = sqrt(mu/ri);
xx0 = [ri; 0; 0; 0; vi; 0; 1000]; % 7x1 vector (also mass included)
fprintf('\nInitial state in inertial J2000 frame [km], [km/s]: \n')
fprintf('%.6f\n', xx0(1:3))
fprintf('%.8f\n', xx0(4:6))

% Target state in J2000 frame
rf = hf + Re;
vf = sqrt(mu/rf);
xxf = [rf; 0; 0; 0; vf*cos(delta_i); vf*sin(delta_i)]; % 6x1 vector
fprintf('\nTarget state in inertial J2000 frame [km], [km/s]: \n')
fprintf('%.6f\n', xxf(1:3))
fprintf('%.8f\n', xxf(4:end))


%% EX 3.2: Adimensionalize the problem

% TU and VU imposing mu_bar = 1
VU = vi/sqrt(DU/ri); % [km/s]
TU = DU/VU; % [s]

% Adimensionalized parameters (ad)
ri_ad = ri/DU;
rf_ad = rf/DU;
vi_ad = vi/VU;
vf_ad = vf/VU;
m0_ad = m0/MU;
Tmax_ad = Tmax/(MU*1000*VU/TU);
Is_ad = Is/TU;
g0_ad = g0/(VU/TU);
rho0_ad = rho0/DU;

% Adimensional initial state
xx0_ad = NaN(7,1);
xx0_ad(1:3) = xx0(1:3)./DU;
xx0_ad(4:6) = xx0(4:6)./VU;
xx0_ad(end) = m0_ad;

% Adimensional target state
xxf_ad = NaN(6,1);
xxf_ad(1:3) = xxf(1:3)./DU;
xxf_ad(4:6) = xxf(4:6)./VU;


%% EX 3.3: PMP problem

% Derived by hand and described in the report


%% EX 3.4: Solve the problem

%--------------------------------------------------------------------------
% The while loop to find X0_best (lambda0 + tf) is commented out
% The solution X0_best of the while loop is used afterward 
%--------------------------------------------------------------------------

% % WHILE LOOP
%
% % Lambda0 + tf
% X0 = NaN(8,1); 
% 
% % Random guesses
% X0_best = [500*rand()-250; 500*rand()-250; 500*rand()-250; 500*rand()-250; 500*rand()-250; 500*rand()-250; 250*rand(); 2*pi*(10 + rand() - 0.5)];
% 
% % Initialize while loop
% Nmax = 10;
% iter = 0;
% fval_temporary = 1;
% exitflag = -1;
% perc = 0.1;
% t0 = 0;
% rv_f = xxf_ad;
% 
% % Options for first iteration
% options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter', 'FunctionTolerance',1e-8,'StepTolerance',1e-10, 'MaxFunctionEvaluations',1e4,'MaxIterations',5e3);
% 
% while iter < Nmax && exitflag ~= 1
% 
%     % First iteration
%     if iter == 0  
% 
%         % Guess lambda0 7x1 vector
%         lambda0_guess = X0_best(1:7);
% 
%         % Build the 8x1 initial guess for the fsolve
%         X0(1:7) = lambda0_guess;
%         X0(8) = X0_best(8);
% 
%         % Solve
%         [X0_opt, fval, exitflag] = fsolve(@(X0) zero_find_prob(X0, xx0_ad,t0, rv_f, Tmax_ad, Is_ad, g0_ad, k1, k2, rho0_ad),X0,options);
% 
%         % Update new best initial guess
%         if norm(fval) < fval_temporary 
%             fval_temporary = norm(fval);
%             X0_best = X0_opt;
%             fval_best = fval;
%             exit_best = exitflag;
%         end
%     end
% 
%     % Successive iterations
%     if iter > 0
% 
%         % Guess lambda0 7x1 vector by perturbing an old solution
%         lambda0_guess = X0_best(1:7);
%         for k = 1 : size(lambda0_guess,1)
%             lambda0_guess(k) = lambda0_guess(k) + perc*lambda0_guess(k)*(-1+2*rand());
%         end
% 
%         % Build  the 8x1 initial guess for the fsolve by perturbing tf
%         X0(1:7) = lambda0_guess; 
%         X0(8) = X0_best(8) + perc*X0_best(8)*(-1+2*rand());
% 
%         % Solve
%         [X0_opt, fval, exitflag] = fsolve(@(X0) zero_find_prob(X0, xx0_ad,t0, rv_f, Tmax_ad, Is_ad, g0_ad, k1, k2, rho0_ad),X0,options);
% 
%         % Update new best initial guess
%         if norm(fval) < fval_temporary 
%             fval_temporary = norm(fval);
%             X0_best = X0_opt;
%             fval_best = fval;
%             exit_best = exitflag;
%         end        
%     end
% 
%     iter = iter + 1;
% 
%     % Reset the algorithm to default after the first iteration
%     options = optimoptions('fsolve','Display','iter', 'FunctionTolerance',1e-8,'StepTolerance',1e-10, 'MaxFunctionEvaluations',1e4,'MaxIterations',5e3);
% end

% Solution of the while loop
X0_best = 1.0e+02 * [-2.149802211903612; -0.103658597224262; 0.008855602315577; 
    -0.103929071057342; -2.146094523457102; -1.129433963233121; 0.025964176697976; 0.644801096581315];
fprintf('\nLambdas: \n')
fprintf('%.4f\n', X0_best(1:7))
fprintf('tf [min]: \n')
fprintf('%.4f\n', X0_best(end)*TU/60)

% Propagation
t0 = 0;
options = odeset('reltol', 3e-14, 'abstol', 3e-14);
X0_opt = NaN(14,1);
X0_opt(1:7) = xx0_ad;
X0_opt(8:14) = X0_best(1:7);
tf = X0_best(8);
[tt, YY] = ode113(@(t,yy) rhs_aug(t, yy, Tmax_ad, Is_ad, g0_ad, k1, k2, rho0_ad), [t0 tf], X0_opt, options);

% Errors evaluation
r_err = (YY(end,1:3)' - xxf_ad(1:3)).*DU;
v_err = (YY(end,4:6)' - xxf_ad(4:6)).*VU;
fprintf('\nPosition error [mm]: \n')
fprintf('%.6f\n', norm(r_err)*1e6)
fprintf('\nVelocity error [micron/s]: \n')
fprintf('%.6f\n', norm(v_err)*1e9)

% Final mass
mf = YY(end,7)*MU;
fprintf('\nFinal mass [kg]: \n')
fprintf('%.4f\n', mf)

% Initial and target orbits
[x_in,y_in] = circular_orbit(ri_ad, 0, 0);
z_in = zeros(size(x_in));
[x_fin,y_fin] = circular_orbit(rf_ad, 0, 0);
z_fin = zeros(size(x_fin));
R_inc = [1, 0, 0;
       0, cos(delta_i), -sin(delta_i);
       0, sin(delta_i),  cos(delta_i)];
rotated_orb = NaN(3,length(z_fin));
for k = 1:length(z_fin)
    rotated_orb(:,k) = R_inc*[x_fin(k) y_fin(k) z_fin(k)]';
end

% 3D plot
figure
plot3(YY(:,1), YY(:,2), YY(:,3),'LineWidth',4);
hold on
plot3(x_in,y_in,z_in,'LineWidth',4, 'Color', "#A2142F")
plot3(rotated_orb(1,:),rotated_orb(2,:),rotated_orb(3,:),'LineWidth',4, 'Color', "#EDB120")
quiver3(0, 0, 0, 1.5, 0, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3); 
quiver3(0, 0, 0, 0, 1.5, 0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
quiver3(0, 0, 0, 0, 0, 0.01, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.3);
grid on
xlabel('x [-]','Interpreter','latex','FontSize',55)
ylabel('y [-]','Interpreter','latex','FontSize',55)
zlabel('z [-]','Interpreter','latex','FontSize',55)
text(1.5,0,0, '$x$', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
text(0,1.5,0, '$y$', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
text(0,0,0.01, '$z$', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',40);
plot3(xx0_ad(1), xx0_ad(2), xx0_ad(3),'diamond','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',20)
plot3(xxf_ad(1), xxf_ad(2), xxf_ad(3),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
legend(["Transfer orbit","Initial orbit","Target orbit","","","","Departure point","Arrival point"],'Location','best','Interpreter','Latex',FontSize=33)
%title('Propagation of optimal solution','Interpreter','Latex',FontSize = 40)

% Polar plot
radius = log(vecnorm(YY(:,1:3)'));
true_anom = atan2(YY(:,2)',YY(:,1)');
figure
polarplot(true_anom,radius,'LineWidth',2)
hold on
polarplot(0,log(ri_ad),'diamond','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',20)
polarplot(0,log(rf_ad),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
%legend(["Transfer orbit","Departure point","Arrival point"],'Location','best','Interpreter','Latex',FontSize=33)
%title('Logaritmic plot of radius module','Interpreter','latex')

% Hamiltonian
H_vect = NaN(size(tt));
for k = 1:length(tt)
    H_vect(k) = hamiltonian(YY(k,:)', Tmax_ad, Is_ad, g0_ad, k1, k2, rho0_ad);
end
figure
plot(tt*TU/60,H_vect,'LineWidth',4)
xlabel('Time [min]','Interpreter','latex','FontSize',55)
ylabel('H(t)','Interpreter','latex','FontSize',55)
xlim([0 tf*TU/60])
grid on
%title('Hamiltonian', 'Interpreter', 'latex');

% Primer vector components in NTW reference frame
lambda_v = YY(:,11:13)';
alpha_ntw = NaN(3,length(tt));
for i = 1:length(tt)
    alpha_inertial = -lambda_v(:,i)/norm(lambda_v(:,i));

    % Rotation matrix from j2000 to NTW frame
    rot_ntw = NaN(3,3);

    T_axis = YY(i,4:6)/norm(YY(i,4:6));
    R_axis = YY(i,1:3)/norm(YY(i,1:3));
    N_axis = cross(R_axis,T_axis);
    W_axis = cross(N_axis,T_axis);

    rot_ntw(:,1) = N_axis;
    rot_ntw(:,2) = T_axis;
    rot_ntw(:,3) = W_axis;

    % Projection of alpha in new frame
    alpha_ntw(:,i) = rot_ntw'*alpha_inertial;
end

% Find points where the declination is zero (i.e AN, DN nodes)
dec = asind(YY(:,3)'./vecnorm(YY(:,1:3)'));
indices = find((dec(1:end-1) < 0 & dec(2:end) > 0) | (dec(1:end-1) > 0 & dec(2:end) < 0));
times_zero_decl = tt(indices)*TU/60;

% Components plot
figure
subplot(3,1,1)
plot(tt.*TU/60,alpha_ntw(1,:),'LineWidth',4)
ylabel('$\alpha^*_N$','Interpreter','latex','FontSize',55)
hold on
for k = 1:length(indices)
    plot([times_zero_decl(k) times_zero_decl(k)],[-1 1],'k--', 'LineWidth', 2)
end
xlim([0 tf*TU/60])
grid on

subplot(3,1,2)
plot(tt.*TU/60,alpha_ntw(2,:),'LineWidth',4)
ylabel('$\alpha^*_T$','Interpreter','latex','FontSize',55)
hold on
for k = 1:length(indices)
    plot([times_zero_decl(k) times_zero_decl(k)],[0 1],'k--', 'LineWidth', 2)
end
xlim([0 tf*TU/60])
grid on

subplot(3,1,3)
plot(tt.*TU/60,alpha_ntw(3,:),'LineWidth',4)
xlabel('Time [min]','Interpreter','latex','FontSize',55)
ylabel('$\alpha^*_W$','Interpreter','latex','FontSize',55)
hold on
for k = 1:length(indices)
    plot([times_zero_decl(k) times_zero_decl(k)],[-0.1 0],'k--', 'LineWidth', 2)
end
xlim([0 tf*TU/60])
grid on

% Check norm is always 1:
check_norm = vecnorm(alpha_ntw);


%% 5) Solve with Tmax = 2.860 N

% Vector of thrust up to Tmax = 2.860 N
Tmax_vect = linspace(Tmax,2.860,10);

% Initialization of numerical continuation
X0_best_new = X0_best;
X0_new = NaN(8,1);
t0 = 0;
rv_f = xxf_ad;
options = optimoptions('fsolve','Display','iter', 'OptimalityTolerance',1e-10, 'FunctionTolerance',1e-8);

% Numerical continuation loop
for j = 1:length(Tmax_vect)
    
    % Update thrust value
    Tmax_ad_new = Tmax_vect(j)/(MU*1000*VU/TU);
    
    % Update initial guess
    X0_new = X0_best_new;
    
    % Solve
    [X0_opt_new, fval_new] = fsolve(@(x) zero_find_prob(x, xx0_ad,t0, rv_f, Tmax_ad_new, Is_ad, g0_ad, k1, k2, rho0_ad),X0_new,options);

    % Update best solution
    X0_best_new = X0_opt_new;
end

% Results 
fprintf('\nLambdas 2.86 N: \n')
fprintf('%.4f\n', X0_opt_new(1:7))
fprintf('\ntf [min] 2.86 N: \n')
fprintf('%.4f\n', X0_opt_new(end)*TU/60)

% Propagation
t0 = 0;
options = odeset('reltol', 3e-14, 'abstol', 3e-14);
Y0_new = NaN(14,1);
Y0_new(1:7) = xx0_ad;
Y0_new(8:14) = X0_opt_new(1:7);
tf_new = X0_opt_new(8);
[tt_new, YY_new] = ode113(@(t,yy) rhs_aug(t, yy, Tmax_ad_new, Is_ad, g0_ad, k1, k2, rho0_ad), [t0 tf_new], Y0_new, options);

% Errors evaluation
r_err_new = (YY_new(end,1:3)' - xxf_ad(1:3)).*DU;
v_err_new = (YY_new(end,4:6)' - xxf_ad(4:6)).*VU;
fprintf('\nPosition error [mm] 2.86 N: \n')
fprintf('%.6f\n', norm(r_err_new)*1e6)
fprintf('\nVelocity error [ micron/s] 2.86 N: \n')
fprintf('%.6f\n', norm(v_err_new)*1e9)

% Final mass
mf_new = YY_new(end,7)*MU;
fprintf('\nFinal mass [kg] 2.86 N: \n')
fprintf('%.4f\n', mf_new)

% Polar plot
radius_new = log(vecnorm(YY_new(:,1:3)'));
true_anom_new = atan2(YY_new(:,2)',YY_new(:,1)');
figure
polarplot(true_anom_new,radius_new,'LineWidth',2)
hold on
polarplot(0,log(ri_ad),'diamond','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',20)
polarplot(0,log(rf_ad),'diamond','MarkerEdgeColor',"#7E2F8E",'MarkerFaceColor',"#7E2F8E",'MarkerSize',20)
%legend(["Transfer orbit","Departure point","Arrival point"],'Location','best','Interpreter','Latex',FontSize=33)
%title('Logaritmic plot of radius module','Interpreter','latex')

% Hamiltonian
H_vect_new = NaN(size(tt_new));
for k = 1:length(tt_new)
    H_vect_new(k) = hamiltonian(YY_new(k,:)', Tmax_ad_new, Is_ad, g0_ad, k1, k2, rho0_ad);
end
figure
plot(tt_new*TU/60,H_vect_new,'LineWidth',4)
xlabel('Time [min]','Interpreter','latex','FontSize',55)
ylabel('H(t)','Interpreter','latex','FontSize',55)
xlim([0 tf_new*TU/60])
grid on
%title('Hamiltonian', 'Interpreter', 'latex');

% Primer vector components in NTW reference frame
lambda_v_new = YY_new(:,11:13)';
alpha_ntw_new = NaN(3,length(tt_new));
for i = 1:length(tt_new)
    alpha_inertial = -lambda_v_new(:,i)/norm(lambda_v_new(:,i));

    % Rotation matrix from j2000 to NTW frame
    rot_ntw = NaN(3,3);

    T_axis = YY_new(i,4:6)/norm(YY_new(i,4:6));
    R_axis = YY_new(i,1:3)/norm(YY_new(i,1:3));
    N_axis = cross(R_axis,T_axis);
    W_axis = cross(N_axis,T_axis);

    rot_ntw(:,1) = N_axis;
    rot_ntw(:,2) = T_axis;
    rot_ntw(:,3) = W_axis;

    % Projection of alpha in new frame
    alpha_ntw_new(:,i) = rot_ntw'*alpha_inertial;
end

% Find points where the declination is zero (i.e AN, DN nodes)
dec_new = asind(YY_new(:,3)'./vecnorm(YY_new(:,1:3)'));
indices_new = find((dec_new(1:end-1) < 0 & dec_new(2:end) > 0) | (dec_new(1:end-1) > 0 & dec_new(2:end) < 0));
times_zero_decl_new = tt_new(indices_new)*TU/60;

% Components plot
figure
subplot(3,1,1)
plot(tt_new.*TU/60,alpha_ntw_new(1,:), 'LineWidth', 2); 
ylabel('$\alpha^*_N$','Interpreter','latex','FontSize',40)
hold on
for k = 1:length(indices_new)
    plot([times_zero_decl_new(k) times_zero_decl_new(k)],[-1 1],'k--', 'LineWidth', 1.2)
end
xlim([0 tf_new*TU/60])
grid on

subplot(3,1,2)
plot(tt_new.*TU/60,alpha_ntw_new(2,:),'LineWidth',2)
ylabel('$\alpha^*_T$','Interpreter','latex','FontSize',40)
hold on
for k = 1:length(indices_new)
    plot([times_zero_decl_new(k) times_zero_decl_new(k)],[0 1],'k--', 'LineWidth', 1.2)
end
xlim([0 tf_new*TU/60])
grid on

subplot(3,1,3)
plot(tt_new.*TU/60,alpha_ntw_new(3,:),'LineWidth',2)
xlabel('Time [min]','Interpreter','latex','FontSize',40)
ylabel('$\alpha^*_W$','Interpreter','latex','FontSize',40)
hold on
for k = 1:length(indices_new)
    plot([times_zero_decl_new(k) times_zero_decl_new(k)],[-0.02 0],'k--', 'LineWidth', 1.2)
end
xlim([0 tf_new*TU/60])
grid on

% Check norm is always 1:
check_norm_new = vecnorm(alpha_ntw_new);


%% Functions

function [x, y] = circular_orbit(r, x_center, y_center)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Generates the coordinates of a circle to represent a circular orbit.
%   The circle is centered at (x_center, y_center) with radius r, and is
%   discretized into 500 points.
%
% Inputs:
%   r        - Orbital radius of the circular orbit [scalar]
%   x_center - X-coordinate of the orbit center [scalar]
%   y_center - Y-coordinate of the orbit center [scalar]
%
% Outputs:
%   x - X-coordinates of the orbit points [500×1]
%   y - Y-coordinates of the orbit points [500×1]
%--------------------------------------------------------------------------

% Define angular samples
theta = linspace(0, 2*pi, 500).';

% Compute Cartesian coordinates
x = x_center + r*cos(theta);
y = y_center + r*sin(theta);

end

function dydt = rhs_aug(~, yy, Tmax_ad, Isp_ad, g0_ad, k1, k2, rho0_ad)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the right-hand side of the augmented system dynamics,
%   including both equations of motion (EoM) and costate dynamics for
%   optimal control of low-thrust orbital transfer with debris avoidance.
%
% Inputs:
%   ~        - Placeholder for time (not explicitly used)
%   yy       - Augmented state vector [14×1]:
%                 [r(3); v(3); m; λ_r(3); λ_v(3); λ_m]
%              where:
%                 r   = position vector (adimensional)
%                 v   = velocity vector (adimensional)
%                 m   = spacecraft mass (adimensional)
%                 λ_r = costates associated with position
%                 λ_v = costates associated with velocity
%                 λ_m = costate associated with mass
%   Tmax_ad  - Adimensional maximum thrust magnitude
%   Isp_ad   - Adimensional specific impulse
%   g0_ad    - Adimensional gravity constant
%   k1, k2   - Debris spatial density parameters
%   rho0_ad  - Reference orbital radius for debris flux (adimensional)
%
% Outputs:
%   dydt - Time derivative of augmented state [14×1]
%             [dr/dt; dv/dt; dm/dt; dλ_r/dt; dλ_v/dt; dλ_m/dt]
%--------------------------------------------------------------------------

% Initialize derivative
dydt = NaN(14,1);

% --- Extract states and costates -----------------------------------------
r  = yy(1:3);     % Position
v  = yy(4:6);     % Velocity
m  = yy(7);       % Mass
lr = yy(8:10);    % Costate of position
lv = yy(11:13);   % Costate of velocity
% lm = yy(14);    % Costate of mass (not directly needed here)

% --- Norms ---------------------------------------------------------------
rnorm  = norm(r);
lvnorm = norm(lv);

% --- Equations of motion -------------------------------------------------
dydt(1:3) = v;                                             % dr/dt
dydt(4:6) = -r./(rnorm^3) - Tmax_ad/m * lv./lvnorm;        % dv/dt
dydt(7)   = -Tmax_ad/(Isp_ad * g0_ad);                     % dm/dt

% --- Costate dynamics ----------------------------------------------------
dydt(8:10)  = -3/(rnorm^5) * dot(r,lv) * r ...             % dλ_r/dt
              + lv./(rnorm^3) ...
              + (2*k1*(rnorm - rho0_ad) * r./rnorm) ./ ( (k2 + (rnorm-rho0_ad)^2)^2 );
dydt(11:13) = -lr;                                         % dλ_v/dt
dydt(14)    = -lvnorm * Tmax_ad / (m^2);                   % dλ_m/dt

end

function H = hamiltonian(yy, Tmax_ad, Isp_ad, g0_ad, k1, k2, rho0_ad)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Evaluates the Hamiltonian for the optimal control problem of
%   low-thrust orbital transfer with debris avoidance. The Hamiltonian
%   includes contributions from the debris flux penalty, the equations of
%   motion, and the mass flow under continuous thrust.
%
% Inputs:
%   yy       - Augmented state vector [14×1]:
%                 [r(3); v(3); m; λ_r(3); λ_v(3); λ_m]
%              where:
%                 r   = position vector (adimensional)
%                 v   = velocity vector (adimensional)
%                 m   = spacecraft mass (adimensional)
%                 λ_r = costates associated with position
%                 λ_v = costates associated with velocity
%                 λ_m = costate associated with mass
%   Tmax_ad  - Adimensional maximum thrust magnitude
%   Isp_ad   - Adimensional specific impulse
%   g0_ad    - Adimensional gravity constant
%   k1, k2   - Debris spatial density parameters
%   rho0_ad  - Reference orbital radius for debris flux (adimensional)
%
% Outputs:
%   H - Value of the Hamiltonian [scalar]
%--------------------------------------------------------------------------

% Extract states and costates
r  = yy(1:3);
v  = yy(4:6);
m  = yy(7);
lr = yy(8:10);
lv = yy(11:13);
lm = yy(14);

% Norms
rnorm  = norm(r);
lvnorm = norm(lv);

% Debris flux penalty term
q = k1 / (k2 + (rnorm - rho0_ad)^2);

% Hamiltonian evaluation
H = q ...                                 % Debris penalty
    + dot(lr, v) ...                      % Costate-position × velocity
    - dot(r, lv) / (rnorm^3) ...          % Gravitational potential term
    - Tmax_ad * lvnorm / m ...            % Thrust contribution
    - lm * Tmax_ad / (Isp_ad * g0_ad);    % Mass depletion contribution

end

function fun = zero_find_prob(lambda0tf, x_i_ad, t0, rv_f, Tmax_ad, Isp_ad, g0_ad, k1, k2, rho0_ad)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Defines the boundary conditions (BCs) and transversality condition for
%   the indirect optimal control problem of a low-thrust transfer with
%   debris avoidance. This function evaluates the nonlinear system of
%   equations whose solution provides the optimal costates and final time.
%
% Inputs:
%   lambda0tf - Initial costates and final time [8×1]:
%                  [λ_r0(3); λ_v0(3); λ_m0; tf]
%   x_i_ad    - Initial adimensional state [7×1]:
%                  [r0(3); v0(3); m0]
%   t0        - Initial time [scalar]
%   rv_f      - Target Cartesian state at tf [6×1]:
%                  [rf(3); vf(3)]
%   Tmax_ad   - Adimensional maximum thrust magnitude
%   Isp_ad    - Adimensional specific impulse
%   g0_ad     - Adimensional gravity constant
%   k1, k2    - Debris spatial density parameters
%   rho0_ad   - Reference orbital radius for debris flux (adimensional)
%
% Outputs:
%   fun - Residual vector [8×1], to be driven to zero by a root-finder:
%            (1:6) Final state error: [rf; vf] - target [6×1]
%            (7)   Final mass costate condition (λ_m(tf) = 0)
%            (8)   Hamiltonian transversality condition (H(tf) = 0)
%--------------------------------------------------------------------------

% Extract final time
tf = lambda0tf(8);

% Build initial augmented state (state + costates)
yy0 = NaN(14,1);
yy0(1:7)   = x_i_ad;          % Initial state
yy0(8:14)  = lambda0tf(1:7);  % Initial costates

% Integration options
options = odeset('reltol', 3e-14, 'abstol', 3e-14);

% Integrate augmented dynamics
[~, yy_prop] = ode113(@(t,yy) rhs_aug(t, yy, Tmax_ad, Isp_ad, g0_ad, k1, k2, rho0_ad), ...
                      [t0 tf], yy0, options);

yy_prop = yy_prop';

% Extract final augmented state
yy_f = yy_prop(:,end);

% Hamiltonian at final time
H_f = hamiltonian(yy_f, Tmax_ad, Isp_ad, g0_ad, k1, k2, rho0_ad);

% Build residual vector
fun = NaN(8,1);
fun(1:6) = yy_f(1:6) - rv_f;   % Final state error
fun(7)   = yy_f(14);           % λ_m(tf) = 0
fun(8)   = H_f;                % Transversality condition
end

