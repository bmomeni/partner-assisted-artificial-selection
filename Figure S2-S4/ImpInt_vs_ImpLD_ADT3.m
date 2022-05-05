% Modeling toxin degradation in a community with an assisting species
% A: assisting population (supports the growth of D)
% D: degrading populations (removes T)
% T: toxin (inhibits A)
% ADT2:     A(c) = A(c-1) + dt*(rA0 - iT*T(c-1)-A(c-1)/KA0)*A(c-1);
clear

%% Initial conditions
T0 = 10; % initial toxin concentration, ug/ml
A0 = 1e5; % initial A density, cells/ml
D0 = 1e5; % initial D density,cells/ml

%% time range
t0 = 0;
tf = 60; % final time, hrs
dt = 0.01;

%% Parameters
rA0 = 0.2; % A growth rate, 1/hr %%set value
KA0 = 1e8; % A carrying capacity %set value
rD0 = 0.22; % max D growth rate, 1/hr %%evolving
KD0 = 3e8; % D carrying capacity %%evolving
iT = 0.003; % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sA = 1e-7; % growth paramter for A in support of D, ml/(cells.hr) %%set value
dD1 = 3e-8; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
dD2 = 3.5e-8; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving

%% Simulating the dynamics
c = 1;
A1(c) = A0;
D1(c) = D0;
T1(c) = T0;
A2(c) = A0;
D2(c) = D0;
T2(c) = T0;
KA1(c) = KA0*(1-iT*T1(c)/rA0);
KD1(c) = KD0*A1(c)/KA0;
KA2(c) = KA0*(1-iT*T2(c)/rA0);
KD2(c) = KD0*A2(c)/KA0;
trng(c) = t0;
t = t0;

while t < tf
    c = c+1;
    t = t + dt;
    trng(c) = t;

    A1(c) = A1(c-1) + dt*(rA0 - iT*T1(c-1))*(1-A1(c-1)/KA1(c-1))*A1(c-1);
    D1(c) = D1(c-1) + dt*min(sA*A1(c-1),rD0)*max(0,(1 - D1(c-1)/KD1(c-1)))*D1(c-1);
    T1(c) = T1(c-1) - dt*dD1*D1(c-1)*T1(c-1);
    KA1(c) = KA0*(1-iT*T1(c)/rA0);
    KD1(c) = KD0*A1(c)/KA0; 
    A2(c) = A2(c-1) + dt*(rA0 - iT*T2(c-1))*(1-A2(c-1)/KA2(c-1))*A2(c-1);
    D2(c) = D2(c-1) + dt*min(sA*A2(c-1),rD0)*max(0,(1 - D2(c-1)/KD2(c-1)))*D2(c-1);
    T2(c) = T2(c-1) - dt*dD2*max(0,(1 - D2(c-1)/KD2(c-1)))*D2(c-1)*T2(c-1);
    KA2(c) = KA0*(1-iT*T2(c)/rA0);
    KD2(c) = KD0*A2(c)/KA0; 
%     if T(c)<1e-3
%         T(c) = 0;
%     end
%     if A(c)<1e-1
%         A(c) = 0;
%     end
%     if D(c)<1e-1
%         D(c) = 0;
%     end
    dt = 0.05/max([abs(dD1*D1(c-1)) abs(rA0 - iT*T1(c-1)) abs(sA*A1(c-1))]);
end

%% Plot results
figure
semilogy(trng,A1,'color',[0.4 0.1 0.4])
hold on
semilogy(trng,D1,'color',[0.8 0.6 0.2])
semilogy(trng,A2,':','color',[0.4 0.1 0.4])
semilogy(trng,D2,':','color',[0.8 0.6 0.2])
legend('A,II','D,II','A,LD','D,LD')
xlabel('Time (hrs)')
ylabel('Density (cells/ml)')
ylim([5e4 5e8])

figure
plot(trng,T1,'color',[0 0.4 0.1])
hold on
plot(trng,T2,':','color',[0 0.4 0.1])
xlabel('Time (hrs)')
ylabel('Toxin conc. (ug/ml)')
ylim([0 11])

% [m, ind] = min(abs(T - T0/2));
% th = interp1(T(ind-10:ind+10),trng(ind-10:ind+10),T0/2);

