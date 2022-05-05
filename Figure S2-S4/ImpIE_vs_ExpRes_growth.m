% Modeling toxin degradation in a community with an assisting species
% A: assisting population (supports the growth of D)
% D: degrading populations (removes T)
% T: toxin (inhibits A)

clear

%% Initial conditions
T0 = 10; % initial toxin concentration, ug/ml
A0 = 1e5; % initial A density, cells/ml
D0 = 1e5; % initial D density,cells/ml
R0 = 0; % initial conc. of assisting resource

%% time range
t0 = 0;
tf = 120; % final time, hrs
dt = 0.001;

%% Parameters
rA0 = 0.2; % A growth rate, 1/hr %%set value
KA0 = 1e8; % A carrying capacity %set value
rD0 = 0.22; % max D growth rate, 1/hr %%evolving
KD0 = 3e8; % D carrying capacity %%evolving
iT = 0.003; % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sA = 1e-7; % growth paramter for A in support of D, ml/(cells.hr) %%set value
dD = 3e-9; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
% gA = 0.46; % rate of R production by A
% cD = 1; % rate of R consumption by D
% KR = 10; % Monod coefficient for growth of D on R
% gA = 1.1; % rate of R production by A
% cD = 1; % rate of R consumption by D
% KR = 100; % Monod coefficient for growth of D on R
gA = 0.2; % rate of R production by A, fmole/cell/hr
cD = 0.07; % rate of R consumption by D, fmole/cell
KR = 2e5; % Monod coefficient for growth of D on R, fmole/ml

%% Simulating the dynamics
c = 1;
A(c) = A0;
D(c) = D0;
R(c) = R0;
De(c) = D0;
KD(c) = KD0; %*A(c)/KA0;
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
    R(c) = R(c-1) + gA*dt*(rA0 - iT*T(c-1)-A(c-1)/KA0)*A(c-1) - dt*cD*rD0*R(c-1)/(R(c-1)+KR)*De(c-1);
    if R(c)<0 % consumtion exceeds production
        R(c) = 0;
        De(c) = De(c-1);% + min(dt*rD0*R(c)/(R(c)+KR)*De(c-1),(R(c-1)/dt + gA*(rA0 - iT*T(c-1)-A(c-1)/KA0)*A(c-1))/(cD*rD0*R(c-1)/(R(c-1)+KR)));
    else % production exceeds consumption
        De(c) = De(c-1) + dt*rD0*R(c)/(R(c)+KR)*De(c-1);
    end
    T(c) = T(c-1) - dt*dD*D(c-1)*T(c-1);
    Te(c) = Te(c-1) - dt*dD*De(c-1)*Te(c-1);
    KD(c) = KD0*A(c)/KA0; 
    if T(c)<1e-3
        T(c) = 0;
    end
    if Te(c)<1e-3
        Te(c) = 0;
    end
    if A(c)<1e-1
        A(c) = 0;
    end
    if D(c)<1e-1
        D(c) = 0;
    end
    dt = 0.001/max([abs(dD*D(c-1)) abs(rA0 - iT*T(c-1)) abs(sA*A(c-1))]);
end

%% Plot results
figure
semilogy(trng,A)
hold on
semilogy(trng,D)
semilogy(trng,De)   
legend('A','D_i_m','D_e_x')
xlabel('Time (hrs)')
ylabel('Population density (cells/ml)')
xlim([0 100])

figure
plot(trng,R,'k:')
hold on
plot(smooth(trng,15),smooth(R,15),'k')
xlabel('Time (hrs)')
ylabel('Resource conc. (ug/ml)')
xlim([0 100])

figure
plot(trng,T,'color',[0 0.4 0.1])
hold on
plot(trng,Te,'--','color',[0 0.4 0.1])
xlabel('Time (hrs)')
ylabel('Toxin conc. (ug/ml)')
ylim([0 11])
xlim([0 100])

% [m, ind] = min(abs(T - T0/2));
% th = interp1(T(ind-10:ind+10),trng(ind-10:ind+10),T0/2);
