% Modeling toxin degradation in a community with an assisting species
% A: assisting population (supports the growth of D)
% D: degrading populations (removes T)
% T: toxin (inhibits A)

% Modified so that only growing D cells contribute to detox; matches ExpEnz
% with large enz decay rate

clear

%% Initial conditions
T0 = 10; % initial toxin concentration, ug/ml
A0 = 1e5; % initial A density, cells/ml
D0 = 1e5; % initial D density,cells/ml
E0 = 0; % initial enzyme conc.

%% time range
t0 = 0;
tf = 96; % final time, hrs
dt = 0.01;

%% Parameters
rA0 = 0.2; % A growth rate, 1/hr %%set value
KA0 = 1e8; % A carrying capacity %set value
rD0 = 0.22; % max D growth rate, 1/hr %%evolving
KD0 = 3e8; % D carrying capacity %%evolving
iT = 0.003; % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sA = 1e-7; % growth paramter for A in support of D, ml/(cells.hr) %%set value
% dD = 5e-9; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
% dD = 0.9e-8; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
% dD = 1.1e-8; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
dD = 5e-9; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
eE = 2.5e-6; % enzyme production rate per cell, uU/ml/cell
% decE = 0.2; % rate of enzyme decay, 1/hr
% decE = 0.06; % rate of enzyme decay, 1/hr
% decE = 0.02; % rate of enzyme decay, 1/hr
decE = 0.5; % rate of enzyme decay, 1/hr
dE = 1e-3; % degradation parameter for E in removal of T, ml/(uU/hr) %%evolving

%% Simulating the dynamics
c = 1;
A(c) = A0;
D(c) = D0;
E(c) = E0;
KD(c) = KD0*A(c)/KA0;
T(c) = T0;
Te(c) = T0;
trng(c) = t0;
t = t0;

while t < tf
    c = c+1;
    t = t + dt;
    trng(c) = t;

    A(c) = A(c-1) + dt*(rA0 - iT*T(c-1)-A(c-1)/KA0)*A(c-1);
    D(c) = D(c-1) + dt*min(sA*A(c-1),rD0)*max(0,(1 - D(c-1)/KD(c-1)))*D(c-1);
    E(c) = E(c-1) + dt*eE*max(0,(1 - D(c-1)/KD(c-1)))*D(c-1) - dt*decE*E(c-1);
    T(c) = T(c-1) - dt*dD*max(0,(1 - D(c-1)/KD(c-1)))*D(c-1)*T(c-1);
    Te(c) = Te(c-1) - dt*dE*E(c-1)*Te(c-1);
    KD(c) = KD0*A(c)/KA0; 
    if T(c)<1e-3
        T(c) = 0;
    end
    if A(c)<1e-1
        A(c) = 0;
    end
    if D(c)<1e-1
        D(c) = 0;
    end
    dt = 0.005/max([abs(dD*D(c-1)) abs(rA0 - iT*T(c-1)) abs(sA*A(c-1))]);
end

%% Plot results
figure
semilogy(trng,A)
hold on
semilogy(trng,D)
legend('A','D')
xlabel('Time (hrs)')
ylabel('Population density (cells/ml)')
xlim([0 100])

figure
plot(trng,T,'color',[0 0.4 0.1])
hold on
plot(trng,Te,':','color',[0 0.4 0.1])
legend('T_I_m_p_L_D','T_E_x_p_E_n_z')
xlabel('Time (hrs)')
ylabel('Toxin conc. (ug/ml)')
xlim([0 100])
ylim([0 10])

figure
plot(trng,E,'m')
xlabel('Time (hrs)')
ylabel('Enzyme conc. (uU/ml)')
xlim([0 100])

% [m, ind] = min(abs(T - T0/2));
% th = interp1(T(ind-10:ind+10),trng(ind-10:ind+10),T0/2);
