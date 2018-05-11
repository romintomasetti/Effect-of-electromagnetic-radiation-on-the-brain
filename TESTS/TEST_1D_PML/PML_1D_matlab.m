function [signal_in_time] = PML

close all;

%% Number of nodes:
NODES = 1000;

%% Spatial:
dx    = 1e-3;

%% Frequency of the signal:
f     = 900e7;

%% Initialization of the domain:
domain = (0:NODES-1).*dx;

%% Initialization of the fields:
Hy = zeros(1, NODES);
Ez = zeros(1, NODES);

%% Physical parameters of free space:
mu0 = 4*pi*1e-7;
eps0 = 8.8541878176e-12;

%% Speed of light:
c_light = sqrt(1/(mu0*eps0));

%% Time step:
dt = dx/c_light;

%% Initialization of parameters of the grid
mu  = mu0 *ones(1, NODES);
eps = eps0*ones(1, NODES);

%% Initialization of the values causing dissipation
sigma_elec = zeros(1, NODES);

%% Order of the PML:
order = 0;

%% Size of the PML in nodes:
PML_width_nodes = 100;

%% Required reflection:
neta      = sqrt(mu0/eps0); 
R         = 1e-8;
sigma_max = -(order+1)*log(R)/(2*neta*PML_width_nodes*dx);

%% Parameter Sigma_m:
SIGMA_M = 20;

%% "Conductivity":
sigma_elec(NODES-PML_width_nodes:NODES) = ...
        SIGMA_M*((1:1:PML_width_nodes+1)./PML_width_nodes).^order;
sigma_mag = sigma_elec .* mu./eps;

%% Parameters of the gaussian pulse:
period = 1 / f * 10;
MEAN = 2*period;
STD  = period / 4;

%% Computation of the coefficients for H:
a = (1 - sigma_mag*dt./(2*mu)) ./ (1 + sigma_mag*dt./(2*mu));
b = (dt./(mu*dx)) ./ (1 + sigma_mag*dt./(2*mu));

%% Computation of the coefficients for E:
c = (1 - sigma_elec*dt./(2*eps)) ./ (1 + sigma_elec*dt./(2*eps));
d = (dt./(eps*dx)) ./ (1 + sigma_elec*dt./(2*eps));

%% Initializing the figure:
figure;
subplot(2,1,1);
hold on;
plot_1 = plot(domain,Ez,'o-');
title('Electric field [V/m]');
xlabel('x');
ylabel('Electric field');
subplot(2,1,2);
hold on;
plot_2 = plot(domain,Hy,'o-');
title('Magnetic field [A/m]');
xlabel('x');
ylabel('Magnetic field');

%% Initialization of the algorithm
current_time = 0;

for t=0:1900
    
    %% Keeping the previous electric field:
    Ez_prev = Ez;
    
    %% Set the source
    Ez(1) = sin(2*pi*f*current_time)*exp(-(current_time-MEAN)^2/(2*STD^2));
    
    %% Update magnetic field:
    Hy(1:NODES-1) = a(1:NODES-1).*Hy(1:NODES-1) + b(1:NODES-1).*(Ez(2:NODES) - Ez(1:NODES-1));

    %% Update electric field:
    Ez(2:NODES-1) = c(2:NODES-1).*Ez(2:NODES-1) + d(2:NODES-1).*(Hy(2:NODES-1) - Hy(1:NODES-2));
    
    %% Apply ABC boundary condition at the end:
    Ez(end) = Ez_prev(end-1)...
                      + ((c_light*dt - dx)/(c_light*dt + dx))...
                         *(Ez(end-1) - Ez_prev(end));

    %% Update figure:
    if mod(t,10) == 0
        subplot(2,1,1);
        set(plot_1,'YData',Ez);
        subplot(2,1,2);
        set(plot_2,'YData',Hy);
        drawnow;
    end

    %% Saving last electric field for post-processing:
    signal_in_time(:) = Ez;
   
    %% Increment simulation time:
    current_time = current_time + dt;
    
end
end
