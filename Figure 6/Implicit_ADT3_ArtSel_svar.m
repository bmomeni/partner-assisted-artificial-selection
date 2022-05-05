% Modeling toxin degradation in a community with an assisting species
% A: assisting population (supports the growth of D)
% D: degrading populations (removes T)
% T: toxin (inhibits A)

clear

%% Initial conditions
T0 = 10; % initial toxin concentration, ug/ml
A0 = 1e5; % initial A density, cells/ml
D0 = 1e5; % initial D density,cells/ml

%% time range
t0 = 0;
tf = 44; % final time, hrs
dt = 0.01;

%% Parameters
Ne = 10000; % ensemble size
rA0m = 0.2; % A growth rate, 1/hr %%set value
KA0m = 1e8; % A carrying capacity %set value
rD0m = 0.22; % max D growth rate, 1/hr %%evolving
KD0m = 3e8; % D carrying capacity %%evolving
iTm = 0.003; % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sAm = 1e-7; % growth paramter for A in support of D, ml/(cells.hr) %%set value
dDm = 1e-8; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving

svar_rng = linspace(0,0.2,21); % variation in parameters; sigma/mu

nsv = 0;
for svar = svar_rng
    nsv = nsv+1;
    disp(nsv)
    rA0 = skewnormal(rA0m,svar*rA0m,-3,[1,Ne]); % A growth rate, 1/hr %%set value
    KA0 = random('Normal',KA0m,svar*KA0m,[1,Ne]); % A carrying capacity %set value
    rD0 = skewnormal(rD0m,svar*rD0m,-3,[1,Ne]); % max D growth rate, 1/hr %%evolving
    KD0 = random('Normal',KD0m,svar*KD0m,[1,Ne]); % D carrying capacity %%evolving
    iT = random('Normal',iTm,svar*iTm,[1,Ne]); % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
    sA = random('Normal',sAm,svar*sAm,[1,Ne]); % growth paramter for A in support of D, ml/(cells.hr) %%set value
    dD = random('Normal',dDm,0.2*dDm,[1,Ne]); %random('Normal',dDm,svar*dDm,[1,Ne]); % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
    
    Af = zeros(1,Ne);
    Df = zeros(1,Ne);
    Tf = zeros(1,Ne);
    
    %% Simulating the dynamics
    c = 1;
    A(c,:) = random('Poisson',A0,[1,Ne]);
    D(c,:) = random('Poisson',D0,[1,Ne]);
    T(c,:) = random('Uniform',0.9*T0,1.1*T0,[1,Ne]);
    Ai = A(c,:);
    Di = D(c,:);
    Ti = T(c,:);
    KA(c,:) = KA0.*(1-iT.*T(c,:)./rA0);
    KD(c,:) = KD0./KA0.*A(c,:);
    trng(c) = t0;
    t = t0;
    
    while t < tf
        c = c+1;
        t = t + dt;
        %     if mod(c,200)==1
        %         disp(t)
        %     end
        trng(c) = t;
        
        A(c,:) = A(c-1,:) + dt*(rA0 - iT.*T(c-1)).*(1 - A(c-1,:)./KA(c-1,:)).*A(c-1,:);
        D(c,:) = D(c-1,:) + dt*min(sA.*A(c-1,:),rD0).*max(0,(1 - D(c-1,:)./KD(c-1,:))).*D(c-1,:);
        T(c,:) = T(c-1,:) - dt*dD.*D(c-1,:).*T(c-1,:);
        KD(c,:) = KD0.*A(c,:)./KA0;
        KA(c,:) = KA0.*(1-iT.*T(c,:)./rA0);
        %         T(T<1e-3) = 0;
        %         A(A<1e-1) = 0;
        %         D(D<1e-1) = 0;
        dt = 0.05/max([abs(dD.*D(c-1,:)), abs(rA0 - iT.*T(c-1,:)), abs(rD0.*max(0,(1-D(c-1,:)./KD(c-1,:))))]);
    end
    
    Af = A(c,:);
    Df = D(c,:);
    Tf = T(c,:);
    
    [sP,I] = sort(Af+Df,'descend');
%     [sPc,Ic] = sort(rand(1,Ne),'descend');
    
    DetoxImprov(nsv) = 1/mean(dD)*mean(dD(I(1:round(10*Ne/100))));
    DetoxImprov_sd(nsv) = 1/mean(dD)*std(dD(I(1:round(10*Ne/100))));
    DetoxImprov_ci = bootci(100,@mean,1/mean(dD)*dD(I(1:round(10*Ne/100))));
    DetoxImprov_cil(nsv) = DetoxImprov_ci(1);
    DetoxImprov_ciu(nsv) = DetoxImprov_ci(2);
    
    pMW(nsv) = ranksum(dD(I(1:round(10*Ne/100))),dDm*ones(1,round(10*Ne/100)));
    
    X = [1/dDm*dD', 1/KA0m*Af'+1/KD0m*Df', 1/T0*Tf', 1/A0*Ai', 1/D0*Di', 1/T0*Ti', 1/rA0m*rA0', 1/rD0m*rD0', 1/KA0m*KA0', 1/KD0m*KD0', 1/iTm*iT', 1/sAm*sA'];
    
    [R,P,RLO,RUP]=corrcoef(X);
    RSel(nsv) = R(1,2);
    RLOSel(nsv) = RLO(1,2);
    RUPSel(nsv) = RUP(1,2);
end
figure
errorbar(svar_rng,DetoxImprov,DetoxImprov_sd)
xlabel('Stochasticity')
ylabel('Detox Improvement')
hold on
plot([0 0.2],[1 1],':','color',[0.4 0.4 0.4])
xlim([0 0.2])
ylim([0.6 1.4])

figure
errorbar(svar_rng,DetoxImprov,DetoxImprov-DetoxImprov_cil,DetoxImprov_ciu-DetoxImprov)
xlabel('Stochasticity')
ylabel('Detox Improvement')
hold on
plot([0 0.2],[1 1],':','color',[0.4 0.4 0.4])
xlim([0 0.2])
ylim([0.95 1.3])

figure
errorbar(svar_rng,RSel,RSel-RLOSel,RUPSel-RSel)
xlabel('Stochasticity')
ylabel('Correlation coefficient')
hold on
plot([0 0.2],[0 0],':','color',[0.4 0.4 0.4])
xlim([0 0.2])
ylim([-0.1 1])

disp(pMW)

save('Implicit_ADT3_ArtSel_svar_tf44_Ne10000_Top10.mat')

