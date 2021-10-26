clear all
% % % % % % % % % % % % parameter estimation % % % % % % % %  % % %
global lb ub
%fmincon
x0 = [0.1 10 2 0.3 1 50];
lb = log10([0.01 0.1 2 0.001 0.001 100]);
ub = log10([1    5   2 10    10    400]);
% fmincon(@tetgfp,x0,[],[],[],[],lb,ub)
%pso
n_param=length(lb);
nsol=10;
Solution=zeros(nsol,n_param+1);
    options=pso;
    options.PopulationSize=35;
    options.Generations=150;%increase later to 100-200
    options.StallGenLimit= 20;
    options.TolFun = 1e-3;
    options.Display='iter';
for j=1:nsol
    [x,fval,~]=pso(@tetgfperr,n_param,[],[],[],[],lb,ub,[],options);
    Solution(j,:)=cat(2,x, fval);
end
Solution=sortrows(Solution,n_param+1);
Solution(nsol+1,:)=cat(2,lb, 10);
Solution(nsol+2,:)=cat(2,ub, 10);


% % % % % % % % % % % % steady state calculations % % % % % % % %  % % %
load('0206')
x = Solution(1,1:6);
N = [10 25 50];
for i = 1:3
n = N(i); %copy number
% parameters
A = 10^x(1); % 0.1; %kprod/kdeg for tetR promoter
B = 10^x(2); %10; %kprod/kdeg for gfp promoter
m = 10^x(3); %2; %cooperativity of TetR
K = 10^x(4); %0.3; % half repression constant for TetR:PLtet-O1
ka = 10^x(5); %1; %dissociation const for aTc:TetR
f = 10^x(6); %1/0.07;

aT = logspace(-1.5,3,30); %amount of ATC added
TetT = n*A; %total Tet
a = 0.5*(-(ka + TetT - aT) + sqrt(TetT^2+(ka+aT).^2+2*(ka-aT)*TetT)); %free atc
Tet = TetT*ka./(a+ka); %free Tet
G(i,:) = n*B+n*B*f*(K^m./(K^m+(Tet).^m)); %GFP
end
run fig1d_data.m
figure(2); loglog(aT, G(1,:),'k'); hold on; loglog(f1d(:,1,1),f1d(:,2,1),'ko','MarkerFaceColor','k','MarkerSize',3.75)
loglog(aT, G(2,:),'color',[1 0.75 0]); loglog(f1d(:,1,2),f1d(:,2,2),'color',[1 0.75 0],'marker','o','linestyle','none','MarkerFaceColor',[1 0.75 0],'MarkerSize',3.75)
loglog(aT, G(3,:),'r'); loglog(f1d(:,1,3),f1d(:,2,3),'r','marker','o','linestyle','none','MarkerFaceColor','r','MarkerSize',3.75)
title(num2str(10.^x,3))
ylabel('GFP'); xlabel('ATC')

% % pars = 'kpTet kpGFP coop Tet:DNADiss ATC:TetDiss fold';
% % parindex = strsplit(pars,' ');
% % sol = Solution;%(1:10,:);
% % figure(7)
% % for i = 1:size(sol,2)-1
% %     if sol(end-1,i) ~= sol(end,i)
% % subplot(2,3,i)
% %     hold on
% %     [n, y] = hist(sol(1:end-2,i),linspace(sol(end-1,i),sol(end,i),10));
% %     h = bar(y,n,'hist');
% %     set(h,'FaceColor',[1 0 0],'facealpha',.6,'edgecolor',[0 0 0])
% %     set(gca,'XLim',sol([end-1 end],i))
% %     xlabel(parindex{i})
% %     end
% % end