% Spacecraft Guidance and Navigation
% Assignment # 1
% Author: Matteo Portantiolo

%% SYMS
% syms x y z vx vy vz mu 
% 
% r1 = sqrt((x+mu).^2+y.^2+z.^2);
% r2 = sqrt((x+mu-1).^2+y.^2+z.^2);
% v_sq = vx.^2+vy.^2+vz.^2;
% 
% OM = 1/2*(x.^2+y.^2)+(1-mu)/r1+mu/r2+1/2*mu*(1-mu); %OM=U
% J = 2*OM-v_sq;
% 
% dUdx = x + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(3/2)) - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(3/2));
% dUdy = y - (mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
% dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
% 
% f = [vx;vy;vz;2*vy+dUdx;-2*vx+dUdy;dUdz];
% JJ = jacobian(f,[x y z vx vy vz]); %Assemble the matrix A(t)=dfdx

%% START
clearvars; close all; clc;  

% Setting plot options
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter','latex')
set(0, 'defaultAxesFontSize', 28 ,'defaultAxesFontSizeMode', 'manual');

%% EX 1.1: Lagrangian points

% Gravitational parameter
mu = 0.012150;

% Bodies position
x_E = -mu; y_E = 0;  % Earth
x_M = 1-mu; y_M = 0; % Moon

% Compute distances from bodies 1 and 2 (z=0)
r1 = @(xx) sqrt((xx(1) + mu).^2 + xx(2).^2);
r2 = @(xx) sqrt((xx(1) + mu - 1).^2 + xx(2).^2);

% Compute derivative of the potential
dUdxx = @(xx) [xx(1) - (1-mu)/r1(xx).^3*(mu+xx(1)) + mu/r2(xx).^3*(1-mu-xx(1));
                xx(2) - (1-mu)/r1(xx).^3*xx(2) - mu/r2(xx).^3*xx(2)]; 

% Initial guesses
x0_guess = [0.5, 0;
        2, 0;
        -2, 0;
        0.5, 0.5;
        0.5, -0.5];

% Solve to find Lagrange points
L_points = NaN(5,2);
options = optimoptions('fsolve','OptimalityTolerance',1e-10,'Display','none');
for ii = 1:5
    L_points(ii,:) = fsolve(dUdxx,x0_guess(ii,:),options);
end
fprintf('\nLagrangian points:\n');
display(L_points);

% Contour plot
x = linspace(-1.5,1.5,1000);
y = linspace(-1.5,1.5,1000);
[X,Y] = meshgrid(x,y);
norm_grad = grad_norm(X,Y,mu);

figure
contour(X,Y,norm_grad, 0:0.01:2);
colormap(parula)
hcb = colorbar('eastoutside', 'Ticks', 0:0.2:2); 
hcb.FontSize = 15; 
set(hcb, 'TickLabelInterpreter', 'latex', 'fontsize', 33);
xlabel('x [-]',FontSize=40)
ylabel('y [-]',FontSize=40)
hold on 
grid on
plot(L_points(1,1),L_points(1,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L_points(2,1),L_points(2,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L_points(3,1),L_points(3,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L_points(4,1),L_points(4,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L_points(5,1),L_points(5,2),'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(x_E,y_E,'o','MarkerSize',40,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_M,y_M,'o','MarkerSize',24,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(L_points(1,1),L_points(1,2), '$L_1$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',33); 
text(L_points(2,1),L_points(2,2), '$L_2$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',33); 
text(L_points(3,1),L_points(3,2), '$L_3$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',33); 
text(L_points(4,1),L_points(4,2), '$L_4$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',33);  
text(L_points(5,1),L_points(5,2), '$L_5$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',33); 
text(x_E,y_E, '$Earth$', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(x_M,y_M, '$Moon$', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 

% Check triangular point (computed by hand)
L4 = [1-mu-1/2; sqrt(1-(1/2)^2); 0];
L5 = [1-mu-1/2; -sqrt(1-(1/2)^2); 0];

% Derivative of the potential
dUdx = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x+mu-1)./abs(x+mu-1).^3;
int1 = -2:0.01:-mu-1e-3;
int2 = -mu+1e-3:0.01:1-mu-1e-3;
int3 = 1-mu+1e-3:0.01:2;

% Plot
figure
plot(int1,dUdx(int1),'LineWidth',2,'Color',"#0072BD")
hold on 
plot(int2,dUdx(int2),'LineWidth',2,'Color',"#0072BD")
plot(int3,dUdx(int3),'LineWidth',2,'Color',"#0072BD")
axis([-2 2 -30 30])
plot([-mu -mu],[-30 30],'k--')
plot([1-mu 1-mu],[-30 30],'k--')
grid on 
plot(L_points(1,1),L_points(1,2),'o','MarkerSize',14,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L_points(2,1),L_points(2,2),'o','MarkerSize',14,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L_points(3,1),L_points(3,2),'o','MarkerSize',14,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L4(1),L4(2),'o','MarkerSize',14,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L5(1),L5(2),'o','MarkerSize',14,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(x_E,y_E,'o','MarkerSize',30,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(x_M,y_M,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(L_points(1,1),L_points(1,2), '$L_1$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(L_points(2,1),L_points(2,2), '$L_2$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(L_points(3,1),L_points(3,2), '$L_3$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(L4(1),L4(2), '$L_4$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(L5(1),L5(2), '$L_5$', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(x_E,y_E, '$Earth$', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(x_M,y_M, '$Moon$', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
xlabel('x [-]',FontSize=40)
ylabel('y [-]',FontSize=40) 
legend('Derivative of potential','','','Asymptotes','','Lagrange points','','','','', '', '','Interpreter','Latex',FontSize = 33)

% Jacobi constant
r1 = @(x,y,z) sqrt((x+mu).^2+y.^2+z.^2);
r2 = @(x,y,z) sqrt((x+mu-1).^2+y.^2+z.^2);
v_sq = @(vx,vy,vz) vx.^2+vy.^2+vz.^2;
OM = @(x,y,z) 1/2*(x.^2+y.^2)+(1-mu)/r1(x,y,z)+mu/r2(x,y,z)+1/2*mu*(1-mu);
J = @(x,y,z,vx,vy,vz) 2*OM(x,y,z)-v_sq(vx,vy,vz); 
C_vect = NaN(5,1);
for i = 1:5
    C_vect(i) = J(L_points(i,1),L_points(i,2),0,0,0,0); % C=J(L)
end


%% EX 1.2: Halo orbit

% Initial state
x0 = 1.068792441776;
y0 = 0;
z0 = 0.071093328515;
vx0 = 0;
vy0 = 0.319422926485;
vz0 = 0;

xx0 = [x0;y0;z0;vx0;vy0;vz0];

% Parameters
mu = 0.012150;  
tf = 2;
C = 3.09;

% Propagation without correction
[~,~,~,xx_F]  = propagate(0,xx0,tf,mu,true);  
[~,~,~,xx_B]  = propagate(0,xx0,-tf,mu,true); 

figure
hold on
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',2)
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',2)
plot3(xx0(1),xx0(2),xx0(3),'o','MarkerSize',20,'MarkerEdgeColor',"#0072BD",'MarkerFaceColor',"#0072BD")
plot3(x_M,0,0,'o','MarkerSize',24,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
xlabel('x [-]',FontSize=40)
ylabel('y [-]',FontSize=40)
zlabel('z [-]',FontSize=40)
grid on
legend(["Forward","Backward","Initial position","Moon"],'Location','best')
title('Halo orbit without correction')

% Pseudo-Newton method
err_yf = 1;
err_vxf = 1;    
err_vzf = 1;
err_C = 1;
Nmax    = 50;   
iter    = 0;   
tol     = 1e-12; 

% Set the update to the first guess 
x0_new = x0;
z0_new = z0;
vy0_new = vy0;
dt_new = 0;
yf_ref = 0;
vxf_ref = 0;
vzf_ref = 0;
C_ref = C;

% Set the while loop
mat = NaN(7,7);
while (abs(err_yf)>tol || abs(err_vxf)>tol || abs(err_vzf)>tol || abs(err_C)>tol) && iter < Nmax

    % Perform propagation up to event time
    [xf,PHI,te]  = propagate(0,[x0_new;y0;z0_new;vx0;vy0_new;vz0],tf+dt_new,mu);

    % Actual Jacobi constant and RHS of system (for augmented matrix)
    C_e = J(x0_new,y0,z0_new,vx0,vy0_new,vz0);
    [dxdt] = xyzCR3BP(te,xf, mu);

    % Compute the deviation in the final state
    err_yf = yf_ref-xf(2);
    err_vxf = vxf_ref-xf(4);
    err_vzf = vzf_ref-xf(6);
    err_C= C_ref-C_e;

    % Augmented matrix and simplification with constraints 
    mat(1:6,1:6) = PHI;
    mat(1:6,7) = dxdt;
    mat(7,1:6) = dJdx([x0_new;y0;z0_new;vx0;vy0_new;vz0],mu);
    mat(7,7) = 0;
    mat([1 3 5],:) = [];
    mat(:,[2 4 6]) = [];

    % Error vector
    xff = [err_yf;err_vxf;err_vzf;err_C];

    % Compute the correction
    delta = mat\xff;
    x0_new = x0_new+delta(1);
    z0_new = z0_new+delta(2);
    vy0_new = vy0_new+delta(3);
    dt_new = dt_new+delta(4);

    % Update iteration counter
    iter = iter+1;
end

% New initial state
xx0_new = [x0_new;y0;z0_new;vx0;vy0_new;vz0];
fprintf('\nInitial state:\n');
fprintf('%d\n',xx0_new);

% Propagation of solution
[~,~,te, xx_F, tt_F] = propagate(0,[x0;y0;z0;vx0;vy0_new;vz0], te,mu,false);
[~,~,~, xx_B, tt_B] = propagate(0,[x0;y0;z0;vx0;vy0_new;vz0], -te,mu,false);

% Twice period
[xf,~,te, xx, ~] = propagate(0,xx0_new,2*te,mu,false);
err_x0 = (xf-xx0_new);
max_mean_error = abs(err_x0)./max(abs(xx0_new),abs(xf)); 
fprintf('\nPosition error [km]: %+.3e %+.3e\nVelocity error [km/s]: %+.3e %+.3e\nRelative error [km]: %.3e %.3e %.3e %.3e\n', err_x0(1:2),err_x0(3:4),max_mean_error);

% 2D plot
figure
hold on
plot(xx(:,1),xx(:,2),'LineWidth',2)
plot(xx0_new(1),xx0_new(2),'o','MarkerSize',20,'MarkerEdgeColor',"#0072BD",'MarkerFaceColor',"#0072BD")
plot(L_points(2,1),L_points(2,2),'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(x_M,y_M,'o','MarkerSize',12,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(L_points(2,1),L_points(2,2), '$L_2$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',33); 
text(x_M,y_M, '$Moon$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex','FontSize',33); 
grid on
xlabel('x [-]',FontSize=40)
ylabel('y [-]',FontSize=40)
title('2D Halo orbit with correction')
legend(["Orbit","Initial position"],'Location','best', 'Interpreter', 'latex','FontSize',20); 

% 3D plot
figure
hold on
plot3(xx(:,1),xx(:,2),xx(:,3),'LineWidth',2)
plot3(xx0_new(1),xx0_new(2),xx0_new(3),'o','MarkerSize',20,'MarkerEdgeColor',"#0072BD",'MarkerFaceColor',"#0072BD")
plot3(L_points(2,1),L_points(2,2),0,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot3(x_M,y_M,0,'o','MarkerSize',24,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(xx0_new(1),xx0_new(2),xx0_new(3), '$ x_0^{corrected}$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(L_points(2,1),L_points(2,2),0, '$L_2$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(x_M,y_M,0, '$Moon$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
xlabel('x [-]',FontSize=40)
ylabel('y [-]',FontSize=40)
zlabel('z [-]',FontSize=40)
grid on
%legend(["Halo orbit","Corrected initial state"],'Location','best', 'Interpreter', 'latex','FontSize',33); 
%title('3D Halo orbit with differential correction')


%% EX 1.3: Families of Halo

% Vector of Jacobi constant
C_vect = 3.09:-0.01:3.04;
x0_fam = NaN(6,length(C_vect));
iter = NaN(length(C_vect),1);
tf = NaN(length(C_vect),1);

% Plot of families
figure
for i = 1:length(C_vect)

    % Jacobi constant 
    C = C_vect(i);

    % Initial state
    [x0_fam(:,i),tf(i),iter(i)] = find_Halo(xx0,C,mu);
    xx0 = x0_fam(:,i);

    % Twice period propagation
    [~,~,tf, xx, ~] = propagate(0,x0_fam(:,i),2*tf(i),mu,false);

    % Plot
    hold on
    plot3(xx0(1),xx0(2),xx0(3),'o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k')
    plot3(xx(:,1),xx(:,2),xx(:,3),'LineWidth',2)

end

xlabel('x [-]',FontSize=40)
ylabel('y [-]',FontSize=40)
zlabel('z [-]',FontSize=40)
grid on
plot3(L_points(2,1),L_points(2,2),0,'o','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot3(x_M,y_M,0,'o','MarkerSize',24,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5])
text(L_points(2,1),L_points(2,2),0, '$L_2$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
text(x_M,y_M,0, '$Moon$', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Interpreter', 'latex','FontSize',33); 
legend(["Different $ x_0^{corrected}$","Orbit C = 3.09","","Orbit C = 3.08","","Orbit C = 3.07","","Orbit C = 3.06","","Orbit C = 3.05","","Orbit C = 3.04"],'Location','best','FontSize',33); 


%% FUNCTIONS

function norm_grad = grad_norm(x, y, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the norm of the gradient of the effective potential in the
%   Circular Restricted Three-Body Problem (CR3BP).
%
% Inputs:
%   x     - x-coordinate in the rotating frame [non-dimensional]
%   y     - y-coordinate in the rotating frame [non-dimensional]
%   mu    - Mass parameter, defined as m2 / (m1 + m2) [-]
%
% Outputs:
%   norm_grad - Norm of the gradient of the potential [non-dimensional]
%--------------------------------------------------------------------------

% Distances from the two primary bodies
r1 = sqrt((x + mu).^2 + y.^2);
r2 = sqrt((x + mu - 1).^2 + y.^2);

% Partial derivatives of the effective potential
dUdx = x - (1 - mu)./r1.^3 .* (mu + x) + mu./r2.^3 .* (1 - mu - x);
dUdy = y - (1 - mu)./r1.^3 .* y - mu./r2.^3 .* y;

% Norm of the gradient
norm_grad = sqrt(dUdx.^2 + dUdy.^2);

end

function dxdt = xyzCR3BP_STM(~, xx, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the state derivative for the 3D Circular Restricted Three-Body
%   Problem (CR3BP) together with the State Transition Matrix (STM) dynamics.
%   The system includes the equations of motion and the variational equations
%   (linearization of the dynamics).
%
% Inputs:
%   ~   - Placeholder for time (not used).
%   xx  - Extended state vector (42×1):
%           [x; y; z; vx; vy; vz; Phi(:)]
%          where Phi is the 6×6 State Transition Matrix stored column-wise.
%   mu  - Mass parameter, defined as m2 / (m1 + m2) [-].
%
% Outputs:
%   dxdt - Time derivative of the extended state (42×1).
%          Contains [xdot; ydot; zdot; vxdot; vydot; vzdot; Phidot(:)].
%--------------------------------------------------------------------------

% Extract position and velocity components
x  = xx(1);
y  = xx(2);
z  = xx(3);
vx = xx(4);
vy = xx(5);

% Reshape Phi (STM) from vector to 6×6 matrix
Phi = reshape(xx(7:end), 6, 6);

% Distances from the two primaries
r1 = sqrt((x + mu)^2       + y^2 + z^2);
r2 = sqrt((x + mu - 1)^2   + y^2 + z^2);

% Partial derivatives of the potential function
dUdx = x - (1 - mu)/r1^3*(mu + x) + mu/r2^3*(1 - mu - x);
dUdy = y - (1 - mu)/r1^3*y         - mu/r2^3*y;
dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
     - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);

% Jacobian of the equations of motion (6×6 matrix)
dfdx = [ 0, 0, 0, 1, 0, 0
         0, 0, 0, 0, 1, 0
         0, 0, 0, 0, 0, 1
         (mu - 1)/r1^3 - mu/r2^3 - 3*(mu - 1)*(x + mu)^2/r1^5 + 3*mu*(x + mu - 1)^2/r2^5 + 1, ...
         -3*(mu - 1)*y*(x + mu)/r1^5 + 3*mu*y*(x + mu - 1)/r2^5, ...
         -3*(mu - 1)*z*(x + mu)/r1^5 + 3*mu*z*(x + mu - 1)/r2^5, ...
         0, 2, 0
         -3*(mu - 1)*y*(x + mu)/r1^5 + 3*mu*y*(x + mu - 1)/r2^5, ...
         (mu - 1)/r1^3 - mu/r2^3 - 3*(mu - 1)*y^2/r1^5 + 3*mu*y^2/r2^5 + 1, ...
         -3*(mu - 1)*y*z/r1^5 + 3*mu*y*z/r2^5, ...
        -2, 0, 0
         -3*(mu - 1)*z*(x + mu)/r1^5 + 3*mu*z*(x + mu - 1)/r2^5, ...
         -3*(mu - 1)*y*z/r1^5       + 3*mu*y*z/r2^5, ...
         (mu - 1)/r1^3 - mu/r2^3 - 3*(mu - 1)*z^2/r1^5 + 3*mu*z^2/r2^5, ...
         0, 0, 0 ];

% STM derivative
Phidot = dfdx * Phi;

% Assemble right-hand side of extended system
dxdt = zeros(42, 1);

dxdt(1:3)   = xx(4:6);        % Position derivatives
dxdt(4)     = dUdx + 2*vy;    % Acceleration in x
dxdt(5)     = dUdy - 2*vx;    % Acceleration in y
dxdt(6)     = dUdz;           % Acceleration in z
dxdt(7:end) = Phidot(:);      % Flattened STM derivative

end

function dxdt = xyzCR3BP(~, xx, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the state derivative for the 3D Circular Restricted Three-Body
%   Problem (CR3BP) without the State Transition Matrix (STM). The function
%   integrates only the equations of motion in the rotating frame.
%
% Inputs:
%   ~   - Placeholder for time (not used).
%   xx  - State vector (6×1):
%           [x; y; z; vx; vy; vz]
%          where (x, y, z) are positions and (vx, vy, vz) are velocities
%          in the rotating synodic frame.
%   mu  - Mass parameter, defined as m2 / (m1 + m2) [-].
%
% Outputs:
%   dxdt - Time derivative of the state (6×1):
%          [xdot; ydot; zdot; vxdot; vydot; vzdot].
%--------------------------------------------------------------------------

% Extract state components
x  = xx(1);
y  = xx(2);
z  = xx(3);
vx = xx(4);
vy = xx(5);

% Distances from the two primaries
r1 = sqrt((x + mu)^2     + y^2 + z^2);
r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);

% Partial derivatives of the effective potential
dUdx = x - (1 - mu)/r1^3*(mu + x) + mu/r2^3*(1 - mu - x);
dUdy = y - (1 - mu)/r1^3*y        - mu/r2^3*y;
dUdz = (z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
     - (mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);

% Assemble right-hand side
dxdt = zeros(6,1);
dxdt(1:3) = xx(4:6);    % Position derivatives
dxdt(4)   = dUdx + 2*vy; % Acceleration in x
dxdt(5)   = dUdy - 2*vx; % Acceleration in y
dxdt(6)   = dUdz;        % Acceleration in z

end

function [xf, PHIf, tf, xx, tt] = propagate(t0, x0, tf, mu, varargin)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Propagates the 3D Circular Restricted Three-Body Problem (CR3BP) 
%   dynamics together with the State Transition Matrix (STM). Integration
%   is performed using the ODE78 solver with optional event detection for
%   x-z plane crossings.
%
% Inputs:
%   t0       - Initial time [non-dimensional]
%   x0       - Initial state vector (6×1) [x; y; z; vx; vy; vz]
%   tf       - Final time [non-dimensional]
%   mu       - Mass parameter, defined as m2 / (m1 + m2) [-]
%   varargin - Optional:
%                 evtFlag (logical) = true enables x-z plane crossing events
%
% Outputs:
%   xf    - Final state vector (6×1) at tf
%   PHIf  - Final State Transition Matrix (6×6) at tf
%   tf    - Actual final time reached (can differ if events trigger) [nd]
%   xx    - Full integrated state history (N×42 matrix: state + STM)
%   tt    - Time vector corresponding to xx
%--------------------------------------------------------------------------

% Parse optional event flag
if nargin > 4
    evtFlag = varargin{1};
else
    evtFlag = true;
end

% Time of flight
tof = tf - t0;

% Initialize State Transition Matrix at t0
Phi0 = eye(6);

% Extended initial condition: state + STM
x0Phi0 = [x0; Phi0(:)];

% Integration options with event function
options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12, ...
    'Events', @(t, x) xz_plane_crossing(t, x, evtFlag));

% Perform integration using CR3BP dynamics with STM
[tt, xx] = ode78(@(t, x) xyzCR3BP_STM(t, x, mu), [0 tof], x0Phi0, options_STM);

% Extract final state and STM
xf   = xx(end, 1:6)';                  % Final state vector
PHIf = reshape(xx(end, 7:end), 6, 6);  % Final STM
tf   = tt(end);                        % Final integration time

end

function [value, isterminal, direction] = xz_plane_crossing(~, xx, isTerminal)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Event function for ODE integration that detects crossings of the
%   x-z plane in the CR3BP. The event is triggered when y = 0.
%
% Inputs:
%   ~          - Placeholder for time (not used).
%   xx         - State vector (6×1): [x; y; z; vx; vy; vz].
%   isTerminal - Logical flag:
%                  true  -> stop integration at the first crossing
%                  false -> detect crossing but continue integration
%
% Outputs:
%   value      - Event function value (y-coordinate of state).
%   isterminal - Flag passed to ODE solver (same as input isTerminal).
%   direction  - Event direction (0 = all crossings).
%--------------------------------------------------------------------------

% Event condition: crossing the x-z plane (y = 0)
value      = xx(2);

% Termination flag (stop or not)
isterminal = isTerminal;

% Detect crossings in both directions
direction  = 0;

end

function dJdx = dJdx(xx, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes the gradient of the Jacobi constant with respect to the state
%   vector in the 3D Circular Restricted Three-Body Problem (CR3BP).
%
% Inputs:
%   xx - State vector (6×1):
%          [x; y; z; vx; vy; vz]
%        Positions and velocities in the rotating synodic frame.
%   mu - Mass parameter, defined as m2 / (m1 + m2) [-].
%
% Outputs:
%   dJdx - Gradient of the Jacobi constant with respect to the state (1×6).
%          Order corresponds to [∂J/∂x, ∂J/∂y, ∂J/∂z, ∂J/∂vx, ∂J/∂vy, ∂J/∂vz].
%--------------------------------------------------------------------------

% Extract state components
x  = xx(1);
y  = xx(2);
z  = xx(3);
vx = xx(4);
vy = xx(5);
vz = xx(6);

% Gradient of the Jacobi constant
dJdx = [ 2*x + ((2*mu + 2*x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
             - (mu*(2*mu + 2*x - 2))/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
         2*y - (2*mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) ...
             + (2*y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
         (2*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) ...
             - (2*mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
         -2*vx;
         -2*vy;
         -2*vz ];

% Return as row vector (1×6)
dJdx = dJdx';

end

function [xx0_new, tf, iter] = find_Halo(xx0, C, mu)
%--------------------------------------------------------------------------
% Author: Matteo Portantiolo
% Description:
%   Computes corrected initial conditions for a Halo orbit in the CR3BP
%   using a pseudo-Newton method. The algorithm iteratively adjusts the
%   initial state until symmetry conditions and the target Jacobi constant
%   are satisfied within tolerance.
%
% Inputs:
%   xx0 - Initial guess for the state vector (6×1):
%           [x0; y0; z0; vx0; vy0; vz0]
%   C   - Target Jacobi constant [-]
%   mu  - Mass parameter, defined as m2 / (m1 + m2) [-]
%
% Outputs:
%   xx0_new - Corrected initial state vector (6×1) for Halo orbit
%   tf      - Final time after convergence [non-dimensional]
%   iter    - Number of iterations performed
%--------------------------------------------------------------------------

% Extract initial guess components
x0  = xx0(1);
y0  = xx0(2);
z0  = xx0(3);
vx0 = xx0(4);
vy0 = xx0(5);
vz0 = xx0(6);

% Initial integration time guess
tf = 2;

% Define helper functions for CR3BP energy terms
r1   = @(x,y,z) sqrt((x + mu).^2 + y.^2 + z.^2);
r2   = @(x,y,z) sqrt((x + mu - 1).^2 + y.^2 + z.^2);
v_sq = @(vx,vy,vz) vx.^2 + vy.^2 + vz.^2;
OM   = @(x,y,z) 0.5*(x.^2 + y.^2) + (1 - mu)./r1(x,y,z) + mu./r2(x,y,z) ...
                 + 0.5*mu*(1 - mu);
J    = @(x,y,z,vx,vy,vz) 2*OM(x,y,z) - v_sq(vx,vy,vz);

% Pseudo-Newton iteration setup
err_vxf = 1;     % x-velocity error at final state
err_vzf = 1;     % z-velocity error at final state
err_C   = 1;     % Jacobi constant error
Nmax    = 50;    % Maximum number of iterations
iter    = 0;     % Iteration counter
tol     = 1e-12; % Convergence tolerance

% Initialize correction variables
x0_new  = x0;
z0_new  = z0;
vy0_new = vy0;
dt_new  = 0;
vxf_ref = 0;     % Target vx at symmetry point
vzf_ref = 0;     % Target vz at symmetry point
C_ref   = C;     % Target Jacobi constant

mat = NaN(7,7);  % Augmented system matrix (for Newton correction)

% Iterative correction loop
while (abs(err_vxf) > tol || abs(err_vzf) > tol || abs(err_C) > tol) && iter < Nmax
    
    % Propagate trajectory until symmetry event
    [xf, PHI, te] = propagate(0, [x0_new; y0; z0_new; vx0; vy0_new; vz0], tf + dt_new, mu);
    
    % Current Jacobi constant and dynamics at final state
    C_e   = J(x0_new, y0, z0_new, vx0, vy0_new, vz0);
    dxdt  = xyzCR3BP(te, xf, mu);
    
    % Errors with respect to symmetry conditions and Jacobi constant
    err_vxf = vxf_ref - xf(4);
    err_vzf = vzf_ref - xf(6);
    err_C   = C_ref - C_e;
    
    % Build augmented correction matrix
    mat(1:6,1:6) = PHI;
    mat(1:6,7)   = dxdt;
    mat(7,1:6)   = dJdx([x0_new; y0; z0_new; vx0; vy0_new; vz0], mu);
    mat(7,7)     = 0;
    
    % Remove rows/columns for constrained variables
    mat([1 3 5],:) = [];
    mat(:,[2 4 6]) = [];
    
    % Build error vector
    xff = [-xf(2); -xf(4); -xf(6); err_C];
    
    % Solve for correction
    delta = mat \ xff;
    
    % Apply corrections
    x0_new  = x0_new  + delta(1);
    z0_new  = z0_new  + delta(2);
    vy0_new = vy0_new + delta(3);
    dt_new  = dt_new  + delta(4);
    
    % Increase iteration counter
    iter = iter + 1;
end

% Return corrected initial state
xx0_new = [x0_new; y0; z0_new; vx0; vy0_new; vz0];

end

