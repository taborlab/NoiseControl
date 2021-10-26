% deterministic solution for atc-tetr stochastic model
N = 10;
ktettr = 1e-4;
ktlnR = 6e-2;
kb1 = 2.2e-4;
kd1 = 5e-6;
kf1 = 5.8e-4;
kr1 = 2e-5;
kdil = 3.1e-4;
kmd = 2e-3;
kdif = 6.2e-4;
coeff = 0.028;

ATCext = 0; % input

ksyn = N*ktettr*ktlnR/(kmd) % TetR dimer synthesis rate
kin = coeff*ATCext % proportional to ATC input
Ka = kdil/kb1 % kd1 << kdil approximation applied
KD1 = (kr1+kdil)/kf1; % promoter binding
%% find (a/Ka) as a function of atc-ext & TetR dimer synthesis rate
% solve this inline for a/Ka given atc-ext and ksyn (i.e. N) ==>
T2 = (ksyn/kdil)/(1+a/Ka)
kin - kdif*a-(kdil*a*T2/(a+Ka))*(1+2*a/Ka)  % ends up a cubic


p0 = N/(1 + (T2/(KD1+(T2/KD1)*(kdil/kf1)))*(1+T2/KD1)) % ~ = N/(1+T2/KD1)
p1 = p0*T2/(KD1 + (T2*kdil/(kf1*kd1))) % ~ = p0*(T2/KD1)/(1+ T2/KD1)
p2 = N - (p0+p1); % ~ = p1*T2/KD1

T2_T = ksyn/kdil + p1 + 2*p2
sqrt(T2_T/KD1) % equivalent to TetR_T/Ktet w deterministic [compare ksyn corresponding to N = 10]