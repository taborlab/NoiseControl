% Deterministic model simulation with stochastic equivalent parameters:

% %% semi-analytical
% % parameters used for initial submission
% Ka = 408.8; % AHL-LuxR dissociation constant; molecules
% f = 550; % Maximum fold induction of Plux promoter; fold change
% K = 2.93; % Half maximal activation concentration of LuxR-AHL (Ra); molecules
% Rtmin = 14.5*8; % LuxR0, varies with RBS strength; minimum LuxR0 required for activation of positive feedback ~ K*2/(f+1)
% n = 2; % hill coefficient
% scl = 1;

% % finding parameters closer to those used in stochastic model
phi = 0.6/(3.1e-4);
SF = 1/15;
Ka = 0.5/SF;
f = 916/2;
K = 300/(phi*SF*sqrt(2));
Rtmin = 25/(phi*SF);
n = 2;
scl = 1;
%% simple ODE
% very simple ODE as an alternative to solving simultaneous equations
% more accurate in bistable parameter ranges; faster.

kdil = 3.1e-4;
A = 0;
kr = 0.001; kf1 = (kr)/Ka;

pars = struct('K',K','f',f,'alp',Rtmin*kdil, 'kf1',kf1, 'kr',kr,'At',A);
ahlrange = logspace(-2,2,100);
for i = 1:length(ahlrange)
    pars.At = ahlrange(i)*scl;
    [t,R] = ode15s(@luxsimpleode,[0 50]*3600, [Rtmin 0],{},pars);
    ss(i,:) = R(end,:);
end
figure(8); 
col = [1 0.75 0
    0 0 0 
    0 0.5 1
    0 0.5 0];
loglog(ahlrange, ss(:,1)/ss(1,1),'color',[0 0 1],'linewidth',1.5,'linestyle','-'); hold on;
xlabel('AHL (nM)'); ylabel('GFP (fold)')
% xlim([1e-2 1e2])
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',13, ...
% 'LineWidth', 1,'layer','top');

% plot solutions (normalized) against data
col = [1 0.75 0
    0 0 0 
    0 0.5 1
    0 0.5 0];
run luxrdata % data file
figure(8); 
loglog(B34(:,1)*1000,B34(:,2)/B34(1,2),'o','markersize',6,'linestyle','none','color',col(1,:));hold on
loglog(B64(:,1)*1000,B64(:,2)/B64(1,2),'o','markersize',6,'linestyle','none','color',col(2,:)); 
loglog(B32(:,1)*1000,B32(:,2)/B32(1,2),'o','markersize',6,'linestyle','none','color',col(3,:))
loglog(B31(:,1)*1000,B31(:,2)/B31(1,2),'o','markersize',6,'linestyle','none','color',col(4,:))