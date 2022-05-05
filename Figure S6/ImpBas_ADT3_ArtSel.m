% Modeling toxin degradation in a community with an assisting species
% A: assisting population (supports the growth of D)
% D: degrading populations (removes T)
% T: toxin (inhibits A)

clear

%% Initial conditions
T0 = 10; % initial toxin concentration, ug/ml
A0 = 1e4; % initial A density, cells/ml
D0 = 1e4; % initial D density,cells/ml

%% time range
t0 = 0;
tf = 96; % final time, hrs
dt = 0.01;

%% Parameters
Ne = 10000; % ensemble size
rA0m = 0.2; % A growth rate, 1/hr %%set value
KA0m = 1e8; % A carrying capacity %set value
rD0m = 0.22; % max D growth rate, 1/hr %%evolving
KD0m = 3e8; % D carrying capacity %%evolving
iTm = 0.003; % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sAm = 1e-9; % growth paramter for A in support of D, ml/(cells.hr) %%set value
dDm = 1e-6; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving

svar = 0.02; % variation in parameters; sigma/mu
rA0 = skewnormal(rA0m,svar*rA0m,-3,[1,Ne]); % A growth rate, 1/hr %%set value
KA0 = random('Normal',KA0m,svar*KA0m,[1,Ne]); % A carrying capacity %set value
rD0 = skewnormal(rD0m,svar*rD0m,-3,[1,Ne]); % max D growth rate, 1/hr %%evolving
KD0 = random('Normal',KD0m,svar*KD0m,[1,Ne]); % D carrying capacity %%evolving
iT = random('Normal',iTm,svar*iTm,[1,Ne]); % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sA = random('Normal',sAm,svar*sAm,[1,Ne]); % growth paramter for A in support of D, ml/(cells.hr) %%set value
dD = random('Uniform',0.5*dDm,1.5*dDm,[1,Ne]); %random('Normal',dDm,svar*dDm,[1,Ne]); % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving

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
    
    A(c,:) = A(c-1,:) + dt*(rA0 - iT.*T(c-1)).*(1 - A(c-1,:)./KA(c,:)).*A(c-1,:);
    D(c,:) = D(c-1,:) + dt*sA.*A(c-1,:).*(1 - D(c-1,:)./KD0).*D(c-1,:);
    T(c,:) = T(c-1,:) - dt*dD.*D(c-1,:).*T(c-1,:);
    KA(c,:) = KA0.*(1-iT.*T(c,:)./rA0);
    %         T(T<1e-3) = 0;
    %         A(A<1e-1) = 0;
    %         D(D<1e-1) = 0;
    dt = 0.05/max([abs(dD.*D(c-1,:)), abs(rA0 - iT.*T(c-1,:)), abs(sA.*A(c-1,:).*(1-D(c-1,:)./KD0))]);
end

Af = A(c,:);
Df = D(c,:);
Tf = T(c,:);
%
figure
plot(Af+Df,Tf,'k.')
xlabel('Total cell density')
ylabel('Final toxin conc')
%
% figure
% plot(dD,Tf,'k.')
% xlabel('Detox rate')
% ylabel('Final toxin conc')

figure
plot(Af+Df,dD,'k.')
xlabel('Total cell density')
ylabel('Detox rate')
%
figure
plot(Af+Df,rA0,'k.')
ylabel('A growth rate (1/hr)')
xlabel('Total cell density')

figure
plot(Af+Df,rD0,'k.')
ylabel('D growth rate (1/hr)')
xlabel('Total cell density')

figure
plot(Af+Df,KA0,'k.')
ylabel('A carrying capacity (cells/ml)')
xlabel('Total cell density')

figure
plot(Af+Df,KD0,'k.')
ylabel('D carrying capacity (cells/ml)')
xlabel('Total cell density')

figure
plot(Af+Df,iT,'k.')
xlabel('Total cell density')
ylabel('A inhibition factor')

figure
plot(Af+Df,sA,'k.')
xlabel('Total cell density')
ylabel('A growth factor')

[sP,I] = sort(Af+Df,'descend');

figure
subplot(211)
hist(dD)
xlim([0 2e-6])
xlabel('Detox rate')
ylabel('Count')
title('Initial distribution')
subplot(212)
hist(dD(I(1:round(5*Ne/100))))
xlim([0 2e-6])
xlabel('Detox rate')
ylabel('Count')
title('Distribution after one round of selection')

SelStrength = mean(dD(I(1:round(10*Ne/100))))/mean(dD);

X = [1/dDm*dD', 1/KA0m*Af'+1/KD0m*Df', 1/T0*Tf', 1/A0*Ai', 1/D0*Di', 1/T0*Ti', 1/rA0m*rA0', 1/rD0m*rD0', 1/KA0m*KA0', 1/KD0m*KD0', 1/iTm*iT', 1/sAm*sA'];

[R,P,RLO,RUP]=corrcoef(X);
figure
errorbar(1:12,R(1,:),R(1,:)-RLO(1,:),RUP(1,:)-R(1,:))
set(gca,'XTick',1:12,'XTickLabel',{'d_D','Final A+D','Final T','Initial A','Initial D','Initial T','r_A','r_D','K_A','K_D','i_T','s_A'})
ylabel('Correlation coefficient')
hold on
plot([0 13],[0 0],':','color',[0.4 0.4 0.4])
xlim([0 13])
ylim([-1 1])

figure
errorbar(1:12,R(2,:),R(2,:)-RLO(2,:),RUP(2,:)-R(2,:))
set(gca,'XTick',1:12,'XTickLabel',{'d_D','Final A+D','Final T','Initial A','Initial D','Initial T','r_A','r_D','K_A','K_D','i_T','s_A'})
ylabel('Correlation coefficient')
hold on
plot([0 13],[0 0],':','color',[0.4 0.4 0.4])
xlim([0 13])
ylim([-1 1])

save('ImpBas_ADT4_ArtSel_Ne10000_Top10.mat')

