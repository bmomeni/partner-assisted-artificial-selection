% Modeling toxin degradation in a community with an assisting species
% A: assisting population (supports the growth of D)
% D: degrading populations (removes T)
% T: toxin (inhibits A)

clear

rndseed0 = 3725;
rng(rndseed0,'twister');

%% Initial conditions
T0 = 10; % initial toxin concentration, ug/ml
A0 = 1e5; % initial A density, cells/ml
D0 = 1e5; % initial D density,cells/ml

%% time range
t0 = 0;
NCycles = 10; % number of cycles of selection
% tfrng = [48 48 48 46 44 48]; % final time, hrs
tfrng = 48*ones(1,NCycles); % final time, hrs
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

svar = 0.02; % variation in parameters; sigma/mu
rA0 = skewnormal(rA0m,svar*rA0m,-3,[1,Ne]); % A growth rate, 1/hr %%set value
KA0 = random('Normal',KA0m,svar*KA0m,[1,Ne]); % A carrying capacity %set value
rD0 = skewnormal(rD0m,svar*rD0m,-3,[1,Ne]); % max D growth rate, 1/hr %%evolving
KD0 = random('Normal',KD0m,svar*KD0m,[1,Ne]); % D carrying capacity %%evolving
iT = random('Normal',iTm,svar*iTm,[1,Ne]); % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sA = random('Normal',sAm,svar*sAm,[1,Ne]); % growth paramter for A in support of D, ml/(cells.hr) %%set value
dD = random('Uniform',0.5*dDm,1.5*dDm,[1,Ne]); %random('Normal',dDm,svar*dDm,[1,Ne]); % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving

dD0 = dD; % distribution at the beginning of the selection
fs = 0.1; % fraction of population selected at the end of each cycle
hf = 0.7; % heritability factor for traits

% compiled values over cycles
dDc = dD(1:round(fs*Ne));
rD0c = rD0(1:round(fs*Ne));
KD0c = KD0(1:round(fs*Ne));
cycle = zeros(1,round(fs*Ne)); % cycle label

rA0c = [];
KA0c = [];
iTc = [];
sAc = [];

for nc = 1:NCycles
    
    Af = zeros(1,Ne);
    Df = zeros(1,Ne);
    Tf = zeros(1,Ne);
    
    %% Simulating the dynamics
    c = 1;
    A(c,:) = random('Poisson',A0,[1,Ne]);
    D(c,:) = random('Poisson',D0,[1,Ne]);
    T(c,:) = random('Uniform',0.95*T0,1.05*T0,[1,Ne]); %T0*ones(1,Ne);
    Ai = A(c,:);
    Di = D(c,:);
    Ti = T(c,:);
    KA(c,:) = KA0.*(1-iT.*T(c,:)./rA0);
    KD(c,:) = KD0./KA0.*A(c,:);
    trng(c) = t0;
    t = t0;
    
    while t < tfrng(nc)
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
    DetoxImprov(nc) = mean(dD(I(1:round(fs*Ne))))/mean(dD0);
%     figure
%     plot(Af+Df,dD,'k.')
%     hold on
%     plot(Af(I(1:round(fs*Ne)))+Df(I(1:round(fs*Ne))),dD(I(1:round(fs*Ne))),'r.')
%     figure
%     bar([1 1.3 2 2.3 3 3.3], [mean(Tf(dD<1e-8)),mean(Tf(dD>=1e-8)),mean(rD0(dD<1e-8)),mean(rD0(dD>=1e-8)),mean(1e-8*KD0(dD<1e-8)),mean(1e-8*KD0(dD>=1e-8))])
%     pause
    
    rA0c = [rA0c, rA0(I(1:round(fs*Ne)))];
    KA0c = [KA0c, KA0(I(1:round(fs*Ne)))];
    rD0c = [rD0c, rD0(I(1:round(fs*Ne)))];
    KD0c = [KD0c, KD0(I(1:round(fs*Ne)))];
    iTc = [iTc, iT(I(1:round(fs*Ne)))];
    sAc = [sAc, sA(I(1:round(fs*Ne)))];
    dDc = [dDc, dD(I(1:round(fs*Ne)))];
    cycle = [cycle, nc*ones(1,round(fs*Ne))];
    
    idx = randperm(round(fs*Ne));
    for ii = 1:round(Ne/round(fs*Ne))-1
        idx = [idx, randperm(round(fs*Ne))];
    end
    
    rA0 = skewnormal(rA0m,svar*rA0m,-3,[1,Ne]); % A growth rate, 1/hr %%set value
    KA0 = random('Normal',KA0m,svar*KA0m,[1,Ne]); % A carrying capacity %set value
    rD0 = hf*rD0(I(idx)) + (1-hf)*skewnormal(rD0m,svar*rD0m,-3,[1,Ne]); % max D growth rate, 1/hr %%evolving
    KD0 = hf*KD0(I(idx)) + (1-hf)*random('Normal',KD0m,svar*KD0m,[1,Ne]); % D carrying capacity %%evolving
    iT = random('Normal',iTm,svar*iTm,[1,Ne]); % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
    sA = random('Normal',sAm,svar*sAm,[1,Ne]); % growth paramter for A in support of D, ml/(cells.hr) %%set value
    dD = hf*dD(I(idx)) + (1-hf)*random('Uniform',0.5*dDm,1.5*dDm,[1,Ne]); %random('Normal',dDm,svar*dDm,[1,Ne]); % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
end

figure
plot(0:NCycles,[1 DetoxImprov],'o-')
ylim([0.95 1.3])
xlabel('Number of selection cycles')
ylabel('Detox Improvement')

figure
violinplot(dDc,cycle);
xlabel('Selection cycle')
ylabel('Detox coeff.')

figure
violinplot(rD0c,cycle);
xlabel('Selection cycle')
ylabel('Growth rate of D (1/hr)')

figure
violinplot(KD0c,cycle);
xlabel('Selection cycle')
ylabel('Carrying capacity of D (cells/ml)')

% figure
% violinplot(rA0c,cycle);
% xlabel('Selection cycle')
% ylabel('Growth rate of A (1/hr)')
% 
% figure
% violinplot(KA0c,cycle);
% xlabel('Selection cycle')
% ylabel('Carrying capacity of A (cells/ml)')
% 
% figure
% violinplot(sAc,cycle);
% xlabel('Selection cycle')
% ylabel('Growth support of D by A')
% 
% figure
% violinplot(iTc,cycle);
% xlabel('Selection cycle')
% ylabel('Inhibition of A by T')


save(strcat('Implicit_ADT3_ArtSel_',num2str(NCycles),'Cycles_Top',num2str(fs*100),'_hf',num2str(round(hf*100)),'.mat'))
