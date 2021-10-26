function d = luxsimpleode(t,x,pars)
At = pars.At; f = pars.f; K = pars.K; kdil = 3.1e-4; alp = pars.alp; kf1 = pars.kf1; kr = pars.kr;
Rt = x(1); Ra = x(2);
R = Rt - Ra;
a = At-Ra; n=2;

d(1) = alp*(1+f*(Ra/K)^n)/(1+(Ra/K)^n) - kdil*Rt;
d(2) = kf1*R*a - kr*Ra -kdil*Ra;

d = d';