% module comparison
%% deterministic model
Ka = 408.8; % AHL-LuxR dissociation constant; molecules
f = 550; % Maximum fold induction of Plux promoter; fold change
K = 2.93; % Half maximal activation concentration of LuxR-AHL (Ra); molecules
Rtmin = 14.5; % LuxR0, varies with RBS strength; minimum LuxR0 required for activation of positive feedback ~ K*2/(f+1)
n = 2; % hill coefficient

%plots modules of the model: Ra as a function of Rt (the for loop), Rt as
%a function of Ra (hill function)
%comment out if not required
Rt1 = (10.^(-2:0.1:4))*Rtmin;
Ra2 = 10.^(-4:0.1:4); Rt2 = Rtmin*(1+f*(Ra2/K).^n)./(1+(Ra2/K).^n);

A = 1; % at 1nM external conc
Ra1 = 0.5*((Ka + A + Rt1) - sqrt((Ka + A + Rt1).^2-4*A*Rt1));
figure(7); subplot(2,1,1); loglog(Rt2,Ra2,'r'); hold on; xlabel('R_T'); ylabel('R_a'); title('Modules, at 1 nM AHL')
subplot(2,1,1); loglog(Rt1, Ra1,'b');

%% stochastic model's deterministic module breakup
clear variables
Ka = 0.5;
f = 916;
K = 300;
Rtmin = 25;
n = 2;
phi = 0.6/(3.1e-4);
SF = 1/15;
%plots modules of the model: Ra as a function of Rt (the for loop), Rt as
%a function of Ra (hill function)
%comment out if not required
Rt1 = (10.^(-2:0.1:4))*Rtmin;
Ra2 = 10.^(-4:0.1:4); Rt2 = Rtmin*(1+f*(Ra2/K).^n)./(1+2*(Ra2/K).^n);
figure(7); subplot(2,1,2); loglog(Rt2,Ra2,'r'); hold on; xlabel('R_T'); ylabel('R_a'); title('Modules, at 1 nM AHL')

A = 1; % at 1nM external conc
Ra1 = 0.5*((phi*Ka + SF*A + Rt1) - sqrt((phi*Ka + SF*A + Rt1).^2-4*SF*A*phi*Rt1));
subplot(2,1,2); loglog(Rt1, Ra1,'b');

%% deterministic model with parameters equivalent from stochastic: for the same AHL absolute value (1 nM)
clear variables
phi = 0.6/(3.1e-4);
SF = 1/15;
Ka = 0.5/SF;

f = 916/2;
K = 300/(phi*SF*sqrt(2));
Rtmin = 25/(phi*SF);
n = 2;

%plots modules of the model: Ra as a function of Rt (the for loop), Rt as
%a function of Ra (hill function)
%comment out if not required
Rt1 = (10.^(-2:0.1:4))*Rtmin;
Ra2 = 10.^(-4:0.1:4); Rt2 = Rtmin*(1+(f)*(Ra2/K).^n)./(1+(Ra2/K).^n);
figure(7); subplot(2,1,1); loglog(Rt2,Ra2,'r-.'); hold on; xlabel('R_T'); ylabel('R_a'); title('Modules, at 1 nM AHL')

A = 1; % at 1nM external conc
% Ra1 = 0.5*((phi*Ka + SF*A + Rt1) - sqrt((phi*Ka + SF*A + Rt1).^2-4*SF*A*phi*Rt1));
Ra1 = 0.5*((Ka + A + Rt1) - sqrt((Ka + A + Rt1).^2-4*A*Rt1));
subplot(2,1,1); loglog(Rt1, Ra1,'b-.');
