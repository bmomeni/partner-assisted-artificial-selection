% Modeling toxin degradation in a community with an assisting species
% A: assisting population (supports the growth of D)
% D: degrading populations (removes T)
% T: toxin (inhibits A)

clear

%% Initial conditions
T0 = 10; % initial toxin concentration, ug/ml
NA0 = 15; % number of initial A0 values to screen
ND0 = 15; % number of initial D0 values to screen
A0rng = logspace(1,6,NA0); % initial A density, cells/ml
D0rng = logspace(1,6,ND0); % initial D density,cells/ml

%% time range
t0 = 0;
tf = 60; % final time, hrs
dt = 0.01;

%% Parameters
Ne = 1000; % ensemble size
rA0m = 0.2; % A growth rate, 1/hr %%set value
KA0m = 1e8/5; % A carrying capacity %set value
rD0m = 0.22; % max D growth rate, 1/hr %%evolving
KD0m = 3e8/5; % D carrying capacity %%evolving
iTm = 0.003; % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
sAm = 1e-7; % growth paramter for A in support of D, ml/(cells.hr) %%set value
dDm = 1e-8; % degradation parameter for D in removal of T, ml/(cells/hr) %%evolving

MSS = zeros(NA0,ND0);
SSS = zeros(NA0,ND0);

nA0 = 0;
for A0 = A0rng
    nA0 = nA0 + 1;
    disp(nA0)
    nD0 = 0;
    for D0 = D0rng
        nD0 = nD0 + 1;
        % disp(nD0)
        
        Ni = 10;
        SelStrength = zeros(1,Ni);
        
        for ii = 1:Ni
            svar = 0.2; % variation in parameters; sigma/mu
            rA0 = skewnormal(rA0m,0.1*svar*rA0m,-3,[1,Ne]); % A growth rate, 1/hr %%set value
            KA0 = random('Normal',KA0m,0.1*svar*KA0m,[1,Ne]); % A carrying capacity %set value
            rD0 = skewnormal(rD0m,0.1*svar*rD0m,-3,[1,Ne]); % max D growth rate, 1/hr %%evolving
            KD0 = random('Normal',KD0m,0.1*svar*KD0m,[1,Ne]); % D carrying capacity %%evolving
            iT = random('Normal',iTm,0.1*svar*iTm,[1,Ne]); % inhibitory coefficient of T on A, ml/(ug.hr) %%set value
            sA = sAm; %random('Normal',sAm,0.1*svar*sAm,[1,Ne]); % growth paramter for A in support of D, ml/(cells.hr) %%set value
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
                KA(c,:) = KA0.*(1-iT.*T(c,:)./rA0);
                KD(c,:) = KD0.*A(c,:)./KA0;
%                 T(T<1e-3) = 0;
%                 A(A<1e-1) = 0;
%                 D(D<1e-1) = 0;
                dt = 0.05/max([max(abs(dD.*D(c-1,:))), max(abs(rA0 - iT.*T(c-1,:))), max(abs(rD0.*max(0,(1-D(c-1,:)./KD(c-1,:)))))]);
            end
            
            Af = A(c,:);
            Df = D(c,:);
            Tf = T(c,:);
            %
            % figure
            % plot(Af+Df,Tf,'k.')
            % xlabel('Total cell density')
            % ylabel('Final toxin conc')
            %
            % figure
            % plot(dD,Tf,'k.')
            % xlabel('Detox rate')
            % ylabel('Final toxin conc')
            %
            % figure
            % plot(Af+Df,dD,'k.')
            % xlabel('Total cell density')
            % ylabel('Detox rate')
            
            [sP,I] = sort(Af+Df,'descend');
            
            % figure
            % subplot(211)
            % hist(dD)
            % xlim([0 5e-9])
            % subplot(212)
            % hist(dD(I(1:round(Ne/10))))
            % xlim([0 5e-9])
            nGen(ii) = mean(1/log(2)*log(D(c,:)./D(1,:)));
            SelStrength(ii) = mean(dD(I(1:round(Ne/10))))/mean(dD);
            Tfin(ii) = mean(Tf);
        end
        MNG(nA0,nD0) = mean(nGen);
        MSS(nA0,nD0) = mean(SelStrength);
        SSS(nA0,nD0) = std(SelStrength);
        MTF(nA0,nD0) = mean(Tfin);
    end
end

figure
imagesc(log10(D0rng),log10(A0rng),MSS)
axis equal
axis tight
axis xy
colorbar
xlabel('Initial D density (cells/ml)')
ylabel('Initial A density (cells/ml)')
title('Mean Selection Strength')

figure
imagesc(log10(D0rng),log10(A0rng),MNG)
axis equal
axis tight
axis xy
colorbar
xlabel('Initial D density (cells/ml)')
ylabel('Initial A density (cells/ml)')
title('Mean Number of Generations')

figure
plot(100/T0*(T0-MTF),MSS,'k.')
xlabel('Degradation efficiency (%)')
ylabel('Mean selection strength')

figure
semilogx(100/T0*MTF,MSS,'k.')
xlabel('Mean residual T (%)')
ylabel('Mean detox improvement')
xlim([1e-2 1e2])
set(gca,'XTick',[1e-2 1e-1 1e0 1e1 1e2])

disp(mean(mean(MSS)))

save(strcat('ImplInt_ADT3_ArtSel_ScreenA0D0_tf',num2str(tf),'_N',num2str(Ne),'.mat'))
