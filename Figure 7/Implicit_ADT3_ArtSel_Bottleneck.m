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
Ne = 100; % ensemble size
rA0m = 0.2; % A growth rate, 1/hr %%set value
KA0m = 1e8; % A carrying capacity %set value
rD0m = 0.22; % max D growth rate, 1/hr %%evolving
KD0m = 3e8; % D carrying capacity %%evolving
iTm = 0.003; % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sAm = 1e-7; % growth paramter for A in support of D, ml/(cells.hr) %%set value
dDm = 1e-8; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving

Ni = 100;
fbrng = [1 2 3 5 7.5 10 12.5 15 20 25 30];
Nfb = length(fbrng);
DetoxImprov = zeros(Ni,Nfb);

for jj = 1:Nfb
    fb = fbrng(jj);
    
    for ii = 1:Ni
        svar = 0.2; % variation in parameters; sigma/mu
        rA0 = skewnormal(rA0m,0.1*svar*rA0m,-3,[1,Ne]); % A growth rate, 1/hr %%set value
        KA0 = random('Normal',KA0m,0.1*svar*KA0m,[1,Ne]); % A carrying capacity %set value
        rD0 = skewnormal(rD0m,0.1*svar*rD0m,-3,[1,Ne]); % max D growth rate, 1/hr %%evolving
        KD0 = random('Normal',KD0m,0.1*svar*KD0m,[1,Ne]); % D carrying capacity %%evolving
        iT = random('Normal',iTm,0.1*svar*iTm,[1,Ne]); % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
        sA = random('Normal',sAm,0.1*svar*sAm,[1,Ne]); % growth paramter for A in support of D, ml/(cells.hr) %%set value
        dD = random('Uniform',0.5*dDm,1.5*dDm,[1,Ne]); %random('Normal',dDm,svar*dDm,[1,Ne]); % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving
        
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
            dt = 0.05/max([abs(dD.*D(c-1,:)), abs(rA0 - iT.*T(c-1,:)), abs(rD0.*max(0,(1-D(c-1,:)./KD(c-1,:))))]);
        end
        
        Af = A(c,:);
        Df = D(c,:);
        Tf = T(c,:);
        
        [sP,I] = sort(Af+Df,'descend');
        DetoxImprov(ii,jj) = mean(dD(I(1:round(fb*Ne/100))))/mean(dD);
    end
end

MDI = mean(DetoxImprov,1);
SDI = std(DetoxImprov,[],1);
figure
errorbar(100./fbrng,MDI,SDI,'o')
xlabel('Bottleneck stringency')
ylabel('Detox improvement')
hold on
xx = linspace(0,110,100);
yy = 1+0.13*xx./(xx+5);
plot(xx,yy)
set(gca,'XTick',0:20:100)
xlim([0 105])

figure
plot(random('Uniform',-0.2,0.2,[Ni,1])+fbrng,DetoxImprov,'.')
xlabel('Top total cell density cases selected (%)')
ylabel('Detox improvement')
xlim([0 32])

figure
plot(1./sqrt(fbrng*Ne),SDI,'o')
xlabel('1/\sigma_b_n')
ylabel('St. deviation of detox improv.')
hold on
xx = linspace(0,0.1,100);
yy = 2.7*xx;
plot(xx,yy)
set(gca,'XTick',0:0.025:0.1)

save(strcat('Implicit_ADT3_ArtSel_Bottleneck_tf44_Ne',num2str(Ne),'_Ni',num2str(Ni),'.mat'))

