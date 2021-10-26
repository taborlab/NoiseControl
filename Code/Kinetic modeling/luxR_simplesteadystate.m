%% semi-analytical
% parameters used for initial submission
Ka = 408.8; % AHL-LuxR dissociation constant; molecules
f = 550; % Maximum fold induction of Plux promoter; fold change
K = 2.93; % Half maximal activation concentration of LuxR-AHL (Ra); molecules
Rtmin = 10*0.45*0.01/3.1e-4; % LuxR0, varies with RBS strength; minimum LuxR0 required for activation of positive feedback ~ K*2/(f+1)
n = 2; % hill coefficient

% plots modules of the model: Ra as a function of Rt (the for loop), Rt as
% a function of Ra (hill function)
% comment out if not required
% Rt1 = (10.^(-2:0.1:4))*Rtmin;
% Ra2 = 10.^(-4:0.1:4); Rt2 = Rtmin*(1+f*(Ra2/K).^n)./(1+(Ra2/K).^n);
% figure(7);  loglog(Rt2,Ra2,'r'); hold on; xlabel('R_T'); ylabel('R_a'); title('Modules')
% r = logspace(-2,3.7,10);
% for i = 1:length(r)
%          A = r(i);
%          Ra1 = 0.5*((Ka + A + Rt1) - sqrt((Ka + A + Rt1).^2-4*A*Rt1));
%          loglog(Rt1, Ra1,'b');
% end

% solving equations 18 & 21 from SI text:
ahlrange = logspace(-2,2,50);
for i = 1:length(ahlrange)
    A = ahlrange(i);
err = @(Ra,Rt) (Ra-0.5*((Ka + A + Rt) - sqrt((Ka + A + Rt)^2-4*A*Rt)))^2 + (Rt- Rtmin*(1+f*(Ra/K)^n)/(1+(Ra/K)^n))^2; % error function: perhaps superfluous, SSE for Ra and Rt both.
Gfp(i,:) = fmincon(@(x) err(x(1),x(2)),[A f*Rtmin],[],[],[],[],[0 Rtmin],[min(A,Rtmin*f) Rtmin*f]); % second column of fmincon output is Rt, normalized Rt transfer function is proportional to GFP transfer function
% alternate optimization/simultaneous nonlinear equation solver
% err2 = @(Ra,Rt) [Ra-0.5*((Ka + A + Rt) - sqrt((Ka + A + Rt)^2-4*A*Rt)) (Rt- Rtmin*(1+f*(Ra/K)^2)/(1+(Ra/K)^2))];
% Gfp(i,:) = lsqnonlin(@(x) err2(x(1),x(2)),[A f*Rtmin],[0 Rtmin],[min(A,Rtmin*f) Rtmin*f])
end

% plot solutions (normalized) against data
col = [1 0.75 0
    0 0 0 
    0 0.5 1
    0 0.5 0];
run luxrdata % data file
figure(8); loglog(ahlrange, Gfp(:,2)/Gfp(1,2),'linewidth',1.5,'color',[0 0 1]); hold on; xlabel('AHL'); ylabel('normalized R_T'); title('dose response')
loglog(B34(:,1)*1000,B34(:,2)/B34(1,2),'o','markersize',6,'linestyle','none','color',col(1,:));hold on
loglog(B64(:,1)*1000,B64(:,2)/B64(1,2),'o','markersize',6,'linestyle','none','color',col(2,:)); 
loglog(B32(:,1)*1000,B32(:,2)/B32(1,2),'o','markersize',6,'linestyle','none','color',col(3,:))
loglog(B31(:,1)*1000,B31(:,2)/B31(1,2),'o','markersize',6,'linestyle','none','color',col(4,:))


loglog(B33(:,1)*1000,B33(:,2)/B33(1,2),'o','markersize',1.75)

% Derivatives of modules (log-log scale): Ra = f(Rt) and Rt = g(Ra), to
% obtain log gain curves. Used to graphically visualize the analytical condition of
% ultrasensitivity.

% Rt1 = 1*Rtmin;
% xT = logspace(-2,4,25); n = 2;
% alpha= (Ka+xT+Rt1);
% beta = (alpha - sqrt(alpha.^2-4*xT.*Rt1));
% rart = Rt1.*(1-(alpha-2*xT)./(sqrt(alpha.^2-4*xT.*Rt1)))./beta;
% Ra1 = 0.5*((Ka + xT + Rt1) - sqrt((Ka + xT + Rt1).^2-4*xT*Rt1));
% y = Ra1/K;
% fra = (2*(f-1)*y.^2)./((1+y.^2).*(1+f*y.^2));
% 
% figure(4); loglog(xT/Ka,rart,'b'); hold on; loglog(xT/Ka, fra,'r'); loglog(xT/Ka, ones(1,length(xT)),'k')
% legend('L(R_a,R_T)','L(F,R_a)'); xlabel('a_T/K_a'); ylabel('Log Gain')
% figure(5); loglog(xT/Ka, rart.*fra); hold on; loglog(xT/Ka, ones(1,length(xT)),'k');
% xlabel('a_T/K_a'); ylabel('L(R_a,R_T).L(F,R_a)'); title('product of log gains')
%% simple ODE
% very simple ODE as an alternative to solving simulaneous equations
% more accurate in bistable parameter ranges; faster.
% Ka = 100; f = 420; K = 5; alp = 1*kdil;
% X = [0.900864407874771,2.61584667887231,-0.878643904047058+1,0.282558582817605,-1.59079155635523+1]; scl = 1;% B31
% X = [2.83174017598842,2.61922673582138,0.687894778856033,1.74137911914448]; %B31#2; scaled up, nM units
% Ka = 10^X(1); f = 10^X(2); 
% K =10^X(3); Rtmin = 1*10^X(4); 
kdil = 3.1e-4;
A = 0;
kr = 0.001; kf1 = (kr)/Ka;

pars = struct('K',K','f',f,'alp',Rtmin*kdil, 'kf1',kf1, 'kr',kr,'At',A);
ahlrange = logspace(-2,2,100);
for i = 1:length(ahlrange)
    pars.At = ahlrange(i);
    [t,R] = ode15s(@luxsimpleode,[0 50]*3600, [Rtmin 0],{},pars);
    ss(i,:) = R(end,:);
end
figure(8); 
col = [1 0.75 0
    0 0 0 
    0 0.5 1
    0 0.5 0];
loglog(ahlrange, ss(:,1)/ss(1,1),'color',[0 0 1],'linewidth',1.5,'linestyle','--'); hold on;
% xlabel('AHL (a.u.)'); ylabel('GFP (fold)')
xlim([1e-2 1e2])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',13, ...
'LineWidth', 1,'layer','top');
%% simple ODE model optimization
% fitting 4 parameters for the simple ODE model
%           % Ka   f   K    Rmin
% x0 = log10([100  420 5  25]);
% lb = log10([10   415 1   1]);
% ub = log10([1000 425 100 100]);
% % options2 = optimoptions('fmincon','PlotFcn',@optimplotfval,'Display','iter');
% % Y = fmincon(@luxerr_simplesteadystate,x0,[],[],[],[],lb,ub,[],options2)
% options = optimoptions('particleswarm','SwarmSize',10,'PlotFcn',@pswplotbestf,'Display','iter');
% options.MaxIterations = 150; options.FunctionTolerance = 1e-3;
% for num = 1:5
%    [sol,fval]= particleswarm(@luxerr_simplesteadystate,4,lb,ub,options);
%    X(num,:)  = [sol,fval]
% end
% X = [0.900864407874771,2.61584667887231,-0.878643904047058+1,0.282558582817605,-1.59079155635523+1]% B31
% X = [0.8780    2.6232    0.0857    0.2996   -0.5818];
% % optimizing b32,64,34
% lb = ([0    0   0]);
% ub = ([10  10 10]);
% options2 = optimoptions('particleswarm','SwarmSize',10,'PlotFcn',@pswplotbestf,'Display','iter');
% options2.MaxIterations = 150; options2.FunctionTolerance = 1e-3;
% X2 = particleswarm(@(g) luxerr_2(g,X),3,lb,ub,options2);