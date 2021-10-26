b34 = 250; b64 = 220; b32 = 50; b31 = 35; b33 = 3.4; %b64 = 10^X(2); b32 = 10^X(3); b31 = 10^X(4); b33 = 10^X(5);
Ka = 0.5; 
kdif = 0.6; 
K = 300; 
pcopy = 20;
kpd = 3.1e-4; 
kmd = 4e-3; 
ktlnG = 0.39;
ktr = 42*kpd; 
kon = 1e-3; koff = 1.2; KD = koff/kon; f = KD/((1+KD)/520 - 1); %KD = 520; koff = kon*KD; 
kdis = 0.0016;
kra = kdis/(Ka);
ktlnR = b64*kmd;
%% Rate constants
p.ktlnG = ktlnG;
p.ktlnR = ktlnR;
p.ktr = ktr;                        
p.kra = kra;
p.kdis = kdis;
p.kmd = kmd;
p.kpd = kpd;
p.kdif = kdif;
p.kin = kdif*00;
p.kon = kon;
p.koff = koff;
p.K = K;
p.f = f;
%% Stochastic
%%% Specify reaction network
% % Simulate a two-state model of gene expression
import Gillespie.*
pfun = @luxrpropensities;
                %p0 p0ON m  R  G  a  Ra
stoich_matrix = [0  0    0  0  1  0  0   % G trl
                 0  0    0  1  0  0  0   % R trl
                 0  0    1  0  0  0  0   % trx
                 0  0    0 -1  0 -1  1   % R-a binding
                 0  0    0  1  0  1 -1   % R-a unbinding
                 0  0   -1  0  0  0  0   % m decay
                 0  0    0 -1  0  0  0   % R decay
                 0  0    0  0 -1  0  0   % G decay
                 0  0    0  0  0 -1  0   % a decay + outflux
                 0  0    0  0  0  0 -1   % Ra decay
                 0  0    0  0  0  1  0   % a influx
                -1  1    0  0  0  0  0   % gene on
                 1 -1    0  0  0  0  0]; % gene off

b = [b34 b64 b32 b31];

ahlmat = [0 logspace(log10(3.92e-05),log10(0.021),30);
    0 logspace(log10(5.12e-05),log10(0.057),30);
    0 logspace(log10(5.54e-04),log10(0.053),30);
    0 logspace(log10(5.54e-04),log10(0.053),30)
    0 logspace(log10(0.063),log10(14.46),30)]*1000/15; % for 1030 sims
     col = [1 0.75 0
            1 0 0 
            0 0.5 1
            0 0.5 0
            0 0 0];
meanval = zeros(size(ahlmat,2),4); cvval = zeros(size(ahlmat,2),4); 

pdtr = makedist('lognormal',ktr,0.15*ktr);
pdmd = makedist('lognormal',kmd,0.15*kmd);
pdpd = makedist('lognormal',kpd,0.15*kpd);
pdpcopy = makedist('lognormal',20,0.15*20);

% for RBS = 1:4
RBS = 4;
p.ktlnR = b(RBS)*kmd;
% initial state: 0 inducer steady state
p.kin = 0;
ahlrange = ahlmat(RBS,:);
p.ktr = log(random(pdtr));
p.kmd = log(random(pdmd));
p.kpd = log(random(pdpd));
x0    = zeros(1,7); x0(1) = round(log(random(pdpcopy)));
ss0 = zeros(10,7);
for n = 1:200
[~,x] = directMethod(stoich_matrix, pfun, [0 10*3600], x0, p);
ss0(n,:) = x(end,:);
clearvars x
end
pre = round(mean(ss0,1));
% post ahl sim
for i = 1:length(ahlrange)
    ahl = ahlrange(i);
    p.kin = kdif*ahl;
    ss = zeros(50,length(x0));
    tic
        for n = 1:500
    p.ktr = log(random(pdtr));
    p.kmd = log(random(pdmd));
    p.kpd = log(random(pdpd));
    pre(1) = round(log(random(pdpcopy)));
            [~,x] = directMethod(stoich_matrix, pfun, [0 6*3600], pre, p,[],1e8); % 6 hr simulation
            ss(n,:) = x(end,:);
            clearvars x
        end
    save(['C:\Users\sdr4\OneDrive\Lab\Karl_Tabor\LuxR-PF\Data files\1030_konfra_RBS' num2str(RBS) '_' num2str(i)],'ss')
    toc
    % optional plots of distributions
%     figure; histogram(ss(:,5));
%     meanval(i,RBS) = mean(ss(:,5));
%     cvval(i,RBS) = std(ss(:,5))/mean(ss(:,5));
%     disp([meanval(:,RBS) cvval(:,RBS)])
    clearvars ss
end
% end
figure(5); subplot(1,2,1); loglog(ahlrange(1:3)*15, meanval(:,RBS)/pre(5),'color',col(RBS,:),'linewidth',2,'marker','o'); hold on
subplot(1,2,2); semilogx(ahlrange(1:3)*15, cvval(:,RBS),'color',col(RBS,:),'linewidth',2,'marker','o'); hold on

%% unpack saved distributions and plot statistics (Supplementary Figure 3)
for k = 1:31
ss = []; RBS = 2;
load(['/Users/satyajitrao/OneDrive/Lab/Karl_Tabor/LuxR-PF/Data files/1030_konfra_RBS' num2str(RBS) '_' num2str(k)])
meanval(k,RBS) = mean(ss(:,5));
cvval(k,RBS) = std(ss(:,5))/mean(ss(:,5));
end

ahlmat = [0 logspace(log10(3.92e-05),log10(0.021),30);
    0 logspace(log10(5.12e-05),log10(0.057),30);
    0 logspace(log10(5.54e-04),log10(0.053),30);
    0 logspace(log10(5.54e-04),log10(0.053),30)
    0 logspace(log10(0.063),log10(14.46),30)]*1000;
col = [1 0.75 0
            1 0 0 
            0 0.5 1
            0 0.5 0
            0 0 0];

figure(5); subplot(1,2,1); loglog(ahlmat(RBS,1:31), meanval(:,RBS)/meanval(1,RBS),'color',col(RBS,:),'linewidth',2,'marker','none'); hold on
subplot(1,2,2); semilogx(ahlmat(RBS,1:31), cvval(:,RBS),'color',col(RBS,:),'linewidth',2,'marker','none'); hold on

run luxrdata.m
    figure(5); subplot(1,2,1); hold on; loglog(B31(:,1)*1000, B31(:,2)/B31(1,2),'o','markersize',4,'linestyle','none','color',col(4,:))
    subplot(1,2,1); hold on; loglog(B34(:,1)*1000, B34(:,2)/B34(1,2),'o','markersize',4,'linestyle','none','color',col(1,:))
    subplot(1,2,1); hold on; loglog(B64(:,1)*1000, B64(:,2)/B64(1,2),'o','markersize',4,'linestyle','none','color',col(2,:)); 
    subplot(1,2,1); hold on; loglog(B32(:,1)*1000, B32(:,2)/B32(1,2),'o','markersize',4,'linestyle','none','color',col(3,:))
%     subplot(1,2,1); hold on; loglog(B33(:,1)*1000, B33(:,2)/B33(1,2),'o','markersize',4,'linestyle','none','color',col(5,:))
    subplot(1,2,1); xlabel('AHL (nM)'); ylabel('mean sfGFP (fold change)'); 
%     subplot(1,2,2); semilogx(ahlrange,cvval,'color',col(3,:),'linewidth',2'); hold on;
    subplot(1,2,2); semilogx(B31(:,1)*1000,B31(:,3),'o','markersize',4,'linestyle','none','color',col(4,:))
    subplot(1,2,2); semilogx(B34(:,1)*1000,B34(:,3),'o','markersize',4,'linestyle','none','color',col(1,:))
    subplot(1,2,2); semilogx(B64(:,1)*1000,B64(:,3),'o','markersize',4,'linestyle','none','color',col(2,:))
    subplot(1,2,2); semilogx(B32(:,1)*1000,B32(:,3),'o','markersize',4,'linestyle','none','color',col(3,:))
%     subplot(1,2,2); semilogx(B33(:,1)*1000,B33(:,3),'o','markersize',4,'linestyle','none','color',col(5,:))
    subplot(1,2,2); xlabel('AHL (nM)'); ylabel('Noise');
%% mRNA statistics -- not for paper
for RBS = 1:4
for k = 1:31
ss = []; 
load(['/Users/satyajitrao/OneDrive/Lab/Karl_Tabor/LuxR-PF/Data files/1030_konfra_RBS' num2str(RBS) '_' num2str(k)])
meanval_mrna(k,RBS) = mean(ss(:,3));
fano(k,RBS) = var(ss(:,3))/mean(ss(:,3));
end
col = [1 0.75 0
            1 0 0 
            0 0.5 1
            0 0.5 0];
end
figure(67); subplot(1,2,1); semilogx(meanval_mrna, fano,'marker','^','linestyle','none','markersize',6,'linewidth',1.5); hold on;
legend('strongest','strong','weak','weakest')
xlabel('<n> mean mRNA copies/cell'); ylabel('Fano, \sigma^2/<n>')
subplot(1,2,2); loglog(meanval_mrna, fano./meanval_mrna,'marker','^','linestyle','none','markersize',6,'linewidth',1.5)
xlabel('<n> mean mRNA copies/cell'); ylabel('Noise \sigma^2/<n>^2')
%% bonus analysis, not needed for paper
% % % sensitivity/ hill coefficient of simulated GFP transfer
% functions
% ahlmat= [0 0.2 0.4 0.7 0.8 1   2 5;
%          0 0.2 0.4 0.7 0.8 1   2 5 ;
%          0 0.5 1   1.4 2   2.5 5 10;
%          0 0.5 1   1.4 2   2.5 5 10];
% % ahlmat = repmat([0 logspace(-1,1,15)],4,1);
% opf = @(p,x) (p(1) + p(2)*x.^p(3)./(p(4)^p(3)+x.^p(3)))/p(1);
% for RBS = 1:4
% % 1. fit a hill function
% errorfn = @(p) norm(log10(opf(p,ahlmat(RBS,:))') - log10(meanval(:,RBS)/meanval(1,RBS)),2);
% % a = lsqcurvefit(opf,[1 520 2 1],ahlmat(RBS,:)',meanval(:,RBS)/meanval(1,RBS),[1 100 1 0],[]) %; m(RBS) = a(3);
% a = fmincon(errorfn,[1 500 4.2 5],[],[],[],[],zeros(1,4),[]); m(RBS) = a(3);
% inp = logspace(-1,1,20); figure(8); loglog(inp, opf(a,inp),'color',col(RBS,:)); hold on; loglog(ahlmat(RBS,:)',meanval(:,RBS)/meanval(1,RBS),'linestyle','none','marker','o','color',col(RBS,:))
% end
% % 
% % % % hill coef fitting for data - only need once
% % % % run C:\Users\Satyajit\OneDrive\Lab\Karl_Tabor\LuxR-PF\luxrdata
% % % % b = fmincon(@(p) norm(log10(opf(p,B34(:,1))) - log10(B34(:,2)/B34(1,2)),2)^2,[1 520 2 1],[],[],[],[],zeros(1,4),[]);n(1) = b(3);
% % % % inp = logspace(-4,-1,20); figure(9); loglog(inp, opf(b,inp),'color',col(1,:)); hold on; loglog(B34(:,1),B34(:,2)/B34(1,2),'color',col(1,:),'marker','o','linestyle','none')
% % % % b = fmincon(@(p) norm(log10(opf(p,B64(:,1))) - log10(B64(:,2)/B64(1,2)),2)^2,[1 520 2 1],[],[],[],[],zeros(1,4),[]);n(2) = b(3);
% % % % inp = logspace(-4,-1,20); figure(9); loglog(inp, opf(b,inp),'color',col(2,:)); hold on; loglog(B64(:,1),B64(:,2)/B64(1,2),'color',col(2,:),'marker','o','linestyle','none')
% % % % b = fmincon(@(p) norm(log10(opf(p,B32(:,1))) - log10(B32(:,2)/B32(1,2)),2)^2,[1 520 2 1],[],[],[],[],zeros(1,4),[]);n(3) = b(3);
% % % % inp = logspace(-4,-1,20); figure(9); loglog(inp, opf(b,inp),'color',col(3,:)); hold on; loglog(B32(:,1),B32(:,2)/B32(1,2),'color',col(3,:),'marker','o','linestyle','none')
% % % % b = fmincon(@(p) norm(log10(opf(p,B31(:,1))) - log10(B31(:,2)/B31(1,2)),2)^2,[1 520 2 1],[],[],[],[],zeros(1,4),[]);n(4) = b(3);
% % % % inp = logspace(-4,-1,20); figure(9); loglog(inp, opf(b,inp),'color',col(4,:)); hold on; loglog(B31(:,1),B31(:,2)/B31(1,2),'color',col(4,:),'marker','o','linestyle','none')
% n = [3.98688274177782,5.28798013700158,9.73897115957263,9.71025401003970]; % hill coef fits for data
% % 
% figure; plot(1:4, m,'k-o','LineWidth',2); hold on;plot(1:4, n,'c-o','LineWidth',2); legend('simulations','experiments')
% xlabel('RBS strength (decreasing order)'); ylabel('Hill Coefficient from fit')