clear
infile = 'preliminary data toxicity rery afg2 61119.txt';

fid = fopen(strcat(infile),'r');

dt = 10/60; % time step, in hours
Nr = 354; %max time from text file

nr = 6; % number of rows
nc = 5; % number of columns

OD = zeros(nr,nc,Nr);
OD = zeros(nr,nc,Nr);
Temperature = zeros(1,Nr);
ch = fscanf(fid,'%s',1);
while strcmp(ch,'Temperature')==0
    ch = fscanf(fid,'%s',1);
end
while strcmp(ch,'Time')==0
    ch = fscanf(fid,'%s',1);
end 
while strcmp(ch,'H12')==0
    ch = fscanf(fid,'%s',1);
end 

n = 0;
while n < Nr
    n = n+1;
    disp(n)
    Ts = fscanf(fid,'%s',1);
    Temperature(n) = str2double(fscanf(fid,'%s',1));
    for i = 1:nr
        for j = 1:nc
            OD(i,j,n) = str2double(fscanf(fid,'%s',1));
        end
    end
end

while strcmp(ch,'Time')==0
    ch = fscanf(fid,'%s',1);
end 
while strcmp(ch,'H12')==0
    ch = fscanf(fid,'%s',1);
end 

n = 0;
while n < Nr
    n = n+1;
    disp(n)
    Ts = fscanf(fid,'%s',1);
    Temperature(n) = str2double(fscanf(fid,'%s',1));
    for i = 1:nr
        for j = 1:nc
            FL(i,j,n) = str2double(fscanf(fid,'%s',1));
        end
    end
end

fclose(fid);

% Calculate the growth rates
bg_OD = mean(mean(OD(6,1:2,3:5)));
bg_FL = mean(mean(FL(6,1:2,3:5)));
MaxOD = max(OD(:,:,1:Nr),[],3) - bg_OD;

gr = zeros(nr,nc);
for q = 1:nr
    for v = 1:nc
        %u = higher threshold for OD
        %l = lower threshold for OD where it crosses "5e-3"
        [mm lv(q,v)] = min(abs((shiftdim(OD(q,v,1:Nr)- bg_OD - 5e-2,2)))); %lower range to start analysis
        if max(OD(q,v,1:Nr)- bg_OD)<0.2, %if MAX OD doesn't get to 0.2 then it shifts lower what range to look at
            [mm uv(q,v)] = min(abs((shiftdim(OD(q,v,1:Nr)- bg_OD - 0.1,2))));
        else
            uv(q,v) = find((shiftdim(OD(q,v,1:Nr)- bg_OD - 0.2,1))>0,1);
        end
        %Growth Rate
        estpoly = polyfit(dt*(lv(q,v):uv(q,v)),log(shiftdim(OD(q,v,lv(q,v):uv(q,v))- bg_OD,1)),1);%analyzes data
        gr(q,v) = estpoly(1);
    end
end 
% gr(6,1:2) = 0;

te = dt*(1:Nr);
for ii = 1:nr
    for jj = 1:nc
        Q = shiftdim(OD(ii,jj,1:Nr)-mean(OD(ii,jj,1:5)),1); % replace with growth curve from data
        [FitQ, p, ~, ~] = fit_logistic(te,Q);
        
        GR2(ii,jj) = p(3); % growth rate values
        CC(ii,jj) = max(OD(ii,jj,1:Nr)-mean(OD(ii,jj,1:5))); % carrying capacity values
    end
end

m_MaxOD = mean(MaxOD(1:5,:),2);
sd_MaxOD = std(MaxOD(1:5,:),[],2);
m2_MaxOD = mean(CC(1:5,:),2);
sd2_MaxOD = std(CC(1:5,:),[],2);

m_gr = mean(gr(1:5,:),2);
sd_gr = std(gr(1:5,:),[],2);
m2_gr = mean(GR2(1:5,:),2);
sd2_gr = std(GR2(1:5,:),[],2);


T = [10, 20, 50, 100, 0];
figure
% errorbar(T,m_MaxOD,sd_MaxOD,'bs-')
% hold on
errorbar(T,m2_MaxOD,sd2_MaxOD,'ro')
xlabel('AFG_2 Conc. (\mug/ml)')
ylabel('Max OD')
ylim([0 1.4])

figure
% errorbar(T,m_gr,sd_gr,'bs-')
% hold on
errorbar(T,m2_gr,sd2_gr,'ro')
xlabel('AFG_2 Conc. (\mug/ml)')
ylabel('R-ery growth rate (1/hr)')
ylim([0 0.4])

% figure
% for i = 1:nr
%     for j = 1:nc
%     semilogy(0:dt:(Nr-1)*dt,shiftdim(OD(i,j,1:Nr)-bg_OD,1),'b')
%     hold on
%     xlabel('Time (hrs)')
%     ylabel('OD')
%     end
% end
% 
% figure
% for i = 1:nr
%     for j = 1:nc
%     semilogy(0:dt:(Nr-1)*dt,shiftdim(FL(i,j,1:Nr)-bg_FL,1),'r')
%     hold on
%     xlabel('Time (hrs)')
%     ylabel('FL (a.u.)')
%     end
% end

figure
imagesc(MaxOD)
title('Carrying capacity')
colorbar

figure
imagesc(CC)
title('Carrying capacity')
colorbar

figure
imagesc(gr)
title('Growth rate')
colorbar

figure
imagesc(GR2)
title('Growth rate')
colorbar

