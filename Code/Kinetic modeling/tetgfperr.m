function [err,G] = tetgfperr(x)
% Simple steady state model for tet-gfp synthetic construct in ecoli. x
% supplies 6 parameters of the model. Output supplies error wrt exp data to
% PSO algorithm
% outputs: err = norm(log10(exp)-log10(sim))
% outputs: G = steady state GFP predicted by model for 3 values of copy
% number 4,25,50.
global lb ub
% clear all
f1d = []; G = zeros(3,length(f1d(:,1,1))); nerr = zeros(3,1);
run fig1d_data.m
if min(logical(lb<=x & ub>=x)==1)
%   sc101 p15a cole1
N = [4 25 50];
% opts = optimoptions('lsqcurvefit','TolFun',1e-9);
for i = 1:3
% i = 1;
n = N(i); %copy number
% parameters
A = 10^x(1); % 0.1; %Amount of TetR per copy number
B = 10^x(2); %10; %Basal amount of GFP per copy number
m = 10^x(3); %2; %cooperativity of TetR
K = 10^x(4); %0.3; % half repression constant for PLtet-O1
ka = 10^x(5); %1; %dissociation const for aTc:TetR
f = 10^x(6); %1/0.07;

% aT = logspace(-3,1,30); %amount of ATC added
aT = f1d(:,1,i);
TetT = n*A; %total Tet
a = 0.5*(-(ka + TetT - aT) + sqrt(TetT^2+(ka+aT).^2+2*(ka-aT)*TetT)); %free atc
Tet = TetT*ka./(a+ka); %free Tet
G(i,:) = n*B+n*B*f*(K^m./(K^m+(Tet).^m)); %GFP
% fun = @(x,y) y.^x(1)./(y.^x(1)+x(2)^x(1));
% X(i,:) = lsqcurvefit(fun,[1 0.1],aT, (G(i,:)-G(i,1))/(G(i,end)-G(i,1)),[0 0],[],opts);
nerr(i) = sum((log10(f1d(:,2,i))-log10(G(i,:)')).^2);
end
err = sum(nerr);
else
    err = 1e8;
end