% Modeling toxin degradation in a community with an assisting species
% A: assisting population (supports the growth of D)
% D: degrading populations (removes T)
% T: toxin (inhibits A)

clear

%% Initial conditions
T0 = 10; % initial toxin concentration, ug/ml
NA0 = 15; % number of initial A0 values to screen
ND0 = 15; % number of initial D0 values to screen
A0 = 1e5; % initial A density, cells/ml
D0 = 1e5; % initial D density,cells/ml

%% time range
t0 = 0;
tfrng = 20:2:100; % final time, hrs
dt = 0.01;
Nt = length(tfrng);

%% Parameters
Ne = 1000; % ensemble size
rA0m = 0.2; % A growth rate, 1/hr %%set value
KA0m = 1e8/5; % A carrying capacity %set value
rD0m = 0.22; % max D growth rate, 1/hr %%evolving
KD0m = 3e8/5; % D carrying capacity %%evolving
iTm = 0.003; % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sAm = 1e-7; % growth paramter for A in support of D, ml/(cells.hr) %%set value
dDm = 1e-8; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving

MSS = zeros(Nt);
SSS = zeros(Nt);

Ni = 50;
DetoxImprov = zeros(1,Ni);

for ii = 1:Ni
    svar = 0.2; % variation in parameters; sigma/mu
    rA0 = skewnormal(rA0m,0.1*svar*rA0m,-3,[1,Ne]); % A growth rate, 1/hr %%set value
    KA0 = random('Normal',KA0m,0.1*svar*KA0m,[1,Ne]); % A carrying capacity %set value
    rD0 = skewnormal(rD0m,0.1*svar*rD0m,-3,[1,Ne]); % max D growth rate, 1/hr %%evolving
    KD0 = random('Normal',KD0m,0.1*svar*KD0m,[1,Ne]); % D carrying capacity %%evolving
    iT = random('Normal',iTm,0.1*svar*iTm,[1,Ne]); % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
    sA = random('Normal',sAm,0.1*svar*sAm,[1,Ne]); % growth paramter for A in support of D, ml/(cells.hr) %%set value
    dD = random('Normal',dDm,svar*dDm,[1,Ne]); % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
    
    Af = zeros(1,Ne);
    Df = zeros(1,Ne);
    Tf = zeros(1,Ne);
    
    %% Simulating the dynamics
    c = 1;
    A(c,:) = random('Poisson',A0,[1,Ne]);
    D(c,:) = random('Poisson',D0,[1,Ne]);
    T(c,:) = random('Uniform',0.9*T0,1.1*T0,[1,Ne]);
    KA(c,:) = KA0.*(1-iT.*T(c,:)./rA0);
    KD(c,:) = KD0./KA0.*A(c,:);
    trng(c) = t0;
    t = t0;
    
    ct = 1;
    while t < max(tfrng)
        c = c+1;
        t = t + dt;
        trng(c) = t;
        
        A(c,:) = A(c-1,:) + dt*(rA0 - iT.*T(c-1)).*(1 - A(c-1,:)./KA(c-1,:)).*A(c-1,:);
        D(c,:) = D(c-1,:) + dt*min(sA.*A(c-1,:),rD0).*max(0,(1 - D(c-1,:)./KD(c-1,:))).*D(c-1,:);
        T(c,:) = T(c-1,:) - dt*dD.*D(c-1,:).*T(c-1,:);
        KA(c,:) = KA0.*(1-iT.*T(c,:)./rA0);
        KD(c,:) = KD0.*A(c,:)./KA0;
        dt = 0.05/max([dDm*D(c-1,:), rA0, rD0]);
        if t > tfrng(ct)
            disp(ct)
            Af(ct,:) = A(c,:);
            Df(ct,:) = D(c,:);
            Tf(ct,:) = T(c,:);
            [sP,I] = sort(Af(ct,:)+Df(ct,:),'descend');
            nGen(ct,ii) = mean(1/log(2)*log(D(c,:)./D(1,:)));
            DetoxImprov(ct,ii) = mean(dD(I(1:round(Ne/10))))/mean(dD);
            Tfin(ct,ii) = mean(Tf(ct,:));
            
            ct = ct+1;
        end
        
    end
    
    
end
MNG = mean(nGen,2);
MSS = mean(DetoxImprov,2);
SSS = std(DetoxImprov,[],2);
MTF = mean(Tfin,2);

figure
errorbar(tfrng,MSS,SSS)
xlabel('Time (hrs)')
ylabel('Mean detox improvement')

figure
plot(tfrng,MNG)
xlabel('Time (hrs)')
ylabel('Mean Number of Generations')

figure
plot(100/T0*(T0-Tfin),DetoxImprov,'k.')
xlabel('Degradation efficiency (%)')
ylabel('Detox improvement')

figure
semilogx(100/T0*Tfin,DetoxImprov,'k.')
xlabel('Residual T (%)')
ylabel('Detox improvement')
xlim([1e-2 1e2])
set(gca,'XTick',[1e-2 1e-1 1e0 1e1 1e2])

save(strcat('ImpInt_ADT3_ArtSel_st_N',num2str(Ne),'.mat'))
