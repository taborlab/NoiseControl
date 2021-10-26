%% Stochastic
% load and configure solver
sbioloadproject tetgfp2
    % set parameter values
    par = [0.000217030812599453,488.362233961864,0.000581084555141380,1.99525243849006e-05,0.000107391066158416,0.000309232319824481,117.120122241026,3.98575731547017e-06,0.0648124623589506,0.298580479108527,1.56212061492133];
    k1f = par(1);
    f = par(2);
    k9f = par(3); 
    k9b = par(4);
    ktettr = par(5);
    kin = par(6);
    cfac = par(7);
    kgfptr = 3.1e-6;
    kpd = 3.1e-4;%par(9);
    kmd = 2e-3; %par(10);
    ktettrl = par(9);
    kgfptrl = 0.09;
    f2 = par(11);

% %     k1f = 10*1.2364e-5;
% %     f = 485;
% %     k9f = 1.6e-3; 
% %     k9b = 1e-5; 
% %     ktettr = 3.162e-4; 
% %     kin = 5.3e-4; 
% %     cfac = 117;
% %     kgfptr = 3.1e-6;
% %     kpd = 3.1e-4;
% %     kmd = 2e-3;
% %     ktettrl = 0.025;
% %     kgfptrl = 0.15;
            % tet-atc binding
        m1.Reactions(1).KineticLaw(1).Parameters(1).Value = k1f; %Tet-atc forward rate 1/molecule.s
            % Tet promoter binding, dissociation
        m1.Reactions(2).KineticLaw(1).Parameters(1).Value = k9f; %1/(molecule^2 s)
        m1.Reactions(2).KineticLaw(1).Parameters(2).Value = k9b+kpd; %1/s
        
        m1.Reactions(17).KineticLaw(1).Parameters(1).Value = k9f; %1/(molecule^2 s)
        m1.Reactions(17).KineticLaw(1).Parameters(2).Value = k9b; %1/s
        m1.Reactions(19).KineticLaw(1).Parameters(1).Value = kpd; %1/(molecule^2 s)
            % gfp active transcription
        m1.Reactions(4).KineticLaw(1).Parameters(1).Value = kgfptr; % transcription from p2
        m1.Reactions(3).KineticLaw(1).Parameters(1).Value = kgfptr*f; %1/s transcription from p0
        m1.Reactions(18).KineticLaw(1).Parameters(1).Value = kgfptr*f2; %1/s transcription from p1/pTet2
            % tet transcription
        m1.Reactions(8).KineticLaw(1).Parameters(1).Value = ktettr; %1/s
            % atc deg + outflux
        m1.Reactions(14).KineticLaw(1).Parameters(1).Value = kin+kpd; %1/s
            %mRNA lifetime
        m1.Reactions(6).KineticLaw(1).Parameters(1).Value = kmd; %1/s
        m1.Reactions(9).KineticLaw(1).Parameters(1).Value = kmd; %1/s
            %translation rate
        m1.Reactions(5).KineticLaw(1).Parameters(1).Value = kgfptrl; %1/s
        m1.Reactions(7).KineticLaw(1).Parameters(1).Value = ktettrl; %1/s
            % INPUT: ATC influx
        m1.Reactions(13).KineticLaw(1).Parameters(1).Value = 0; % molecules/s
% configure solver
configsetObj = getconfigset(m1,'active');
configsetObj.SolverType = 'ssa';% 'impltau';%'expltau'; %
configsetObj.StopTime = 3600*35;
configsetObj.SolverOptions.LogDecimation = 100;
configsetObj.RunTimeOptions.StatesToLog = 'all';
% initial conditions & initializing variables
f1d = []; run fig1d_data.m
N = [10 25 50];
numRuns = 100;
% rownum = randperm(49,5);
ip=[sort(f1d(:,1,1))';
    sort(f1d(:,1,2))';
    sort(f1d(:,1,3))'];
clearvars f1d cole1 sc101 p15a
% h = zeros(numRuns,size(ip,2),length(N));
meanval = zeros(size(ip));
cvval = zeros(size(ip));
speciesNames = {'GFP'}; 
% for k = 1:3
k=2;
    m1.Species(11).InitialAmount = N(k); % copy #
    m1.Species(9).InitialAmount = N(k);
       for j = 1:size(ip,2)
            m1.Reactions(13).KineticLaw(1).Parameters(1).Value = ip(k,j)*cfac*kin; % atc flux molecule/s
            tic
            simdata = sbioensemblerun(m1, numRuns); %,configsetObj,'linear');
            toc
            stochss=[];
            for i = 1:numRuns
            tind = find(simdata(i,1).Time >= 15*3600,1,'first');
            stochss = cat(1,stochss, simdata(i,1).Data(tind:end,:));
            end
            fn = ['C:\Users\Satyajit\OneDrive\Lab\Karl_Tabor\03222018_SSA_simulationoutputs\GFP_', num2str(N(k)),'_copies_',num2str(j)];
            save fn stochss -v7.3
            meanval(k,j) = mean(stochss(:,7));
            cvval(k,j) = std(stochss(:,7))/mean(stochss(:,7));
            clearvars stochss
        end
% end
%% For the paper
load('/Users/satyajitrao/OneDrive/Lab/Karl_Tabor/Final Simulation data/Low Noise Model/Stats_0327.mat')
col = [0 0 0;
    1 0.75 0;
    1 0 0]; % color matrix for plotting
f1d = []; run fig1d_data.m
ip=[sort(f1d(:,1,1))';
    sort(f1d(:,1,2))';
    sort(f1d(:,1,3))'];
figure(2); % GFP mean vs ATC
subplot(2,3,4)
for k = 1:3
loglog(ip(k,:), meanval(k,:)/meanval(1,1),'color',col(k,:),'linewidth',2); hold on; %errorbar(atcamt,mean(h,1),std(h,1),'-s')        
end
loglog(f1d(:,1,1),f1d(:,2,1)/f1d(1,2,1),'ko','MarkerSize',4) %/f1d(1,2,1)
loglog(f1d(:,1,2),f1d(:,2,2)/f1d(1,2,1),'color',[1 0.75 0],'marker','o','linestyle','none','MarkerSize',4)
loglog(f1d(:,1,3),f1d(:,2,3)/f1d(1,2,1),'r','marker','o','linestyle','none','MarkerSize',4)
ylabel('Mean sfGFP (fold change)'); xlabel('ATC (ng/ml)')
figure(2); % CV vs mean GFP
subplot(2,3,5)
for k = 1:3
semilogx(meanval(k,:)/meanval(1,1), cvval(k,:),'color',col(k,:),'linewidth',2); hold on
end
[~,idx] = sort(sc101(:,1));[~,idx2] = sort(p15a(:,1));[~,idx3] = sort(cole1(:,1));
semilogx(sc101(idx,2)/sc101(1,2),sc101(idx,3),'color',col(1,:),'linestyle','none','marker','o','markersize',4);
semilogx(p15a(idx2,2)/sc101(1,2),p15a(idx2,3),'color',col(2,:),'linestyle','none','marker','o','markersize',4);
semilogx(cole1(idx3,2)/sc101(1,2),cole1(idx3,3),'color',col(3,:),'linestyle','none','marker','o','markersize',4);
xlabel('mean sfGFP (fold change)'); ylabel('Coefficient of Variation')
% legend('Sim 10','Sim 25','Sim 50','Exp sc101','Exp p15a','Exp cole1');


set(gcf,'paperunits','points')
set(gcf,'position',[0 0 235*2 218*2])
saveas(gcf,'fig1.png')


figure(2); % CV vs ATC
subplot(2,3,6)
for k = 1:3
semilogx(ip(k,:), cvval(k,:),'color',col(k,:),'linewidth',2); hold on
end
semilogx(f1d(:,1,1),sc101(:,3),'color',col(1,:),'linestyle','none','marker','o','markersize',4);
semilogx(f1d(:,1,2),p15a(:,3),'color',col(2,:),'linestyle','none','marker','o','markersize',4);
semilogx(f1d(:,1,3),cole1(:,3),'color',col(3,:),'linestyle','none','marker','o','markersize',4);
xlabel('ATC (ng/ml)'); ylabel('Coefficient of Variation')
% legend('Sim 10','Sim 25','Sim 50','Exp sc101','Exp p15a','Exp cole1');
%%
% some distributions and time courses to track

figure; % time course, i = species, j = numrun
for j = 1:numRuns
    for i = 1:10
% i = 1;
subplot(2,5,i)
plot(simdata(j,1).Time/3600,simdata(j,1).Data(:,i)); title(simdata(1,1).DataNames(i)); hold on
pause
    end
end

figure; %histogram
dataj = stochss(:,6);
mh=histogram(dataj);
mh.BinWidth = 1;
mh.Normalization = 'pdf';
% mh.DisplayStyle = 'stairs'; 
mh.EdgeColor = [0.7 0.7 0];
title(simdata(1,1).DataNames(i));
hold on;
% end
%% Deterministic
sbioloadproject tetgfp2
% % % setting parameter values
%     k1f = 10*1.2364e-5;
%     f = 485;
%     k9f = 1.6e-3; 
%     k9b = 1e-5; 
%     ktettr = 3.162e-4; 
%     kin = 5.3e-4; 
%     cfac = 117;
%     kgfptr = 3.1e-6;
%     kpd = 3.1e-4;
%     kmd = 2e-3;
%     ktettrl = 0.025;
%     kgfptrl = 0.15;
par = [0.000217030812599453,488.362233961864,0.000581084555141380,1.99525243849006e-05,0.000107391066158416,0.000309232319824481,117.120122241026,3.98575731547017e-06,0.0648124623589506,0.298580479108527,1.56212061492133];
    k1f = par(1);
    f = par(2);
    k9f = par(3); 
    k9b = par(4);
    ktettr = par(5);
    kin = par(6);
    cfac = par(7);
    kgfptr = 1.4e-6; % par(8);
    kpd = 3.1e-4;%par(9);
    kmd = 2e-3; %par(10);
    ktettrl = par(9);
    kgfptrl = 0.09; %par(10);
    f2 = par(11);
    
            % tet-atc binding
        m1.Reactions(1).KineticLaw(1).Parameters(1).Value = k1f; %Tet-atc forward rate 1/molecule.s
            % Tet promoter binding, dissociation
        m1.Reactions(2).KineticLaw(1).Parameters(1).Value = k9f; %1/(molecule^2 s)
        m1.Reactions(2).KineticLaw(1).Parameters(2).Value = k9b+kpd; %1/s
        
        m1.Reactions(17).KineticLaw(1).Parameters(1).Value = k9f; %1/(molecule^2 s)
        m1.Reactions(17).KineticLaw(1).Parameters(2).Value = k9b; %1/s
        m1.Reactions(19).KineticLaw(1).Parameters(1).Value = kpd; %1/(molecule^2 s)
            % gfp active transcription
        m1.Reactions(4).KineticLaw(1).Parameters(1).Value = kgfptr; % transcription from p2
        m1.Reactions(3).KineticLaw(1).Parameters(1).Value = kgfptr*f; %1/s transcription from p0
        m1.Reactions(18).KineticLaw(1).Parameters(1).Value = kgfptr*f2; %1/s transcription from p1/pTet2
            % tet transcription
        m1.Reactions(8).KineticLaw(1).Parameters(1).Value = ktettr; %1/s
            % atc deg + outflux
        m1.Reactions(14).KineticLaw(1).Parameters(1).Value = kin+kpd; %1/s
            %mRNA lifetime
        m1.Reactions(6).KineticLaw(1).Parameters(1).Value = kmd; %1/s
        m1.Reactions(9).KineticLaw(1).Parameters(1).Value = kmd; %1/s
            %translation rate
        m1.Reactions(5).KineticLaw(1).Parameters(1).Value = kgfptrl; %1/s
        m1.Reactions(7).KineticLaw(1).Parameters(1).Value = ktettrl; %1/s
            % INPUT: ATC influx
        m1.Reactions(13).KineticLaw(1).Parameters(1).Value = 0; % molecules/s
% configure solver
configsetObj = getconfigset(m1,'active');
configsetObj.SolverType = 'ode15s'; 
configsetObj.StopTime = 3600*20;
N = [10 25 50];
f1d = []; run fig1d_data.m
h_det = zeros(3,length(f1d(:,1,1)));

ip=[sort(f1d(:,1,1))';
    sort(f1d(:,1,2))';
    sort(f1d(:,1,3))']; %*2.14e-3*nav*1e-6;
ip = repmat([10.^(-1:0.1:4)],3,1);
% for k = 1:3
k = 2;
m1.Species(11).InitialAmount = N(k); % copy #
m1.Species(9).InitialAmount = N(k);

% for j = 1:size(ip,2)
j = length(f1d(:,1,1));
m1.Reactions(13).KineticLaw(1).Parameters(1).Value = ip(k,j)*cfac*kin; % atc flux molecule/s
simdata = sbiosimulate(m1);
speciesNames = {'GFP'};
[t, x] = selectbyname(simdata, speciesNames);
h_det(k,j) = x(end);
% end
% end

figure(4);
loglog(ip(1,:)*2.14, h_det(1,:)/h_det(1,1),'k'); hold on; %
loglog(ip(2,:)*2.14, h_det(2,:)/h_det(1,1),'color',[1 0.75 0]);
loglog(ip(3,:)*2.14, h_det(3,:)/h_det(1,1),'r');
loglog(f1d(:,1,1)*2.14,f1d(:,2,1)/f1d(1,2,1),'ko','MarkerFaceColor','k','MarkerSize',3.75) %/f1d(1,2,1)
loglog(f1d(:,1,2)*2.14,f1d(:,2,2)/f1d(1,2,1),'color',[1 0.75 0],'marker','o','linestyle','none','MarkerFaceColor',[1 0.75 0],'MarkerSize',3.75)
loglog(f1d(:,1,3)*2.14,f1d(:,2,3)/f1d(1,2,1),'r','marker','o','linestyle','none','MarkerFaceColor','r','MarkerSize',3.75)
title('ode15s-simbio-0308'); ylabel('GFP (deterministic)'); xlabel('ATC')

%%
% % load('NOATC_p15a_EXPDATA')
% % load('NOATC_sc101_EXPDATA')
% % load('NOATC_cole1_EXPDATA')
% % 
% % figure; subplot(3,1,2); z = histogram(exp_p15a,'facecolor',[1 0.75 0],'facealpha',0.5); z.BinWidth = 10; z.Normalization = 'pdf'; hold on; plot(gampdf(0:1:max(exp_p15a),3.58, 47.14),'Color', [1 0.75 0]); xlim([0 800])
% % legend('p15a Exp','\Gamma(3.58,47.14)')
% % subplot(3,1,1); z2= histogram(exp_sc101(1:end-1),'facecolor','k'); z2.BinWidth = 10; z2.Normalization = 'pdf'; hold on; plot(gampdf(0:1:600,3, 56.8),'Color', 'k'); xlim([0 800])
% % legend('sc101 Exp','\Gamma(3,56.8)')
% % subplot(3,1,3); z3= histogram(exp_cole1,'facecolor','r'); z3.BinWidth = 10;z3.Normalization = 'pdf'; hold on; plot(gampdf(0:1:600,3.62, 50.5),'Color', 'r'); xlim([0 800])
% % legend('ColE1 Exp','\Gamma(3.62,50.5)')

% % load('NOATC_p15a_EXPDATA')
% % exp_p15a = exp_p15a - 145;
% % exp_p15a = sort(exp_p15a);
% % ind = find(exp_p15a >= 0, 1, 'first');
% % figure;
% % subplot(3,1,2); z = histogram(exp_p15a(ind:end));
% % gd = fitdist(exp_p15a(ind:end),'Gamma');
% % z.Normalization = 'pdf';
% % hold on; plot(gampdf(0:1:480,gd.a, gd.b),'Color', 'r')
% % std(exp_p15a(ind:end))/mean(exp_p15a(ind:end))
% % mean(exp_p15a(ind:end))
% % legend('p15a Exp-Auto',['\Gamma(',num2str(round(gd.a,2)),',',num2str(round(gd.b,2)),')'])
% % 
% % load('NOATC_cole1_EXPDATA')
% % exp_cole1 = exp_cole1 - 202;
% % exp_cole1 = sort(exp_cole1);
% % ind = find(exp_cole1 >= 0, 1, 'first');
% % 
% % subplot(3,1,3); z = histogram(exp_cole1(ind:end));
% % gd = fitdist(exp_cole1(ind:end),'Gamma');
% % z.Normalization = 'pdf';
% % hold on; plot(gampdf(0:1:480,gd.a, gd.b),'Color', 'r')
% % std(exp_cole1(ind:end))/mean(exp_cole1(ind:end))
% % legend('cole1 Exp-Auto',['\Gamma(',num2str(round(gd.a,2)),',',num2str(round(gd.b,2)),')'])
% % 
% % load('NOATC_sc101_EXPDATA')
% % exp_sc101 = exp_sc101 - 163;
% % exp_sc101 = sort(exp_sc101);
% % ind = find(exp_sc101 >= 0, 1, 'first');
% % subplot(3,1,1); z = histogram(exp_sc101(ind:end-1));
% % gd = fitdist(exp_sc101(ind:end-1),'Gamma');
% % z.Normalization = 'pdf';
% % hold on; plot(gampdf(0:1:480,gd.a, gd.b),'Color', 'r')
% % std(exp_sc101(ind:end-1))/mean(exp_sc101(ind:end-1))
% % legend('sc101 Exp-Auto',['\Gamma(',num2str(round(gd.a,2)),',',num2str(round(gd.b,2)),')'])