function a = luxrpropensities(y, p)
% Return reaction propensities given current state x
p0 = y(1);
p0ON = y(2);
m = y(3);
R = y(4);
G = y(5);
a = y(6);
Ra = y(7);

a = [p.ktlnG*m;       % G trl
     p.ktlnR*m;       % R trl
     p.ktr*p0ON;
     p.kra*R*a;       % trx
     p.kdis*Ra;       % R-a binding
     p.kmd*m;
     p.kpd*R;
     p.kpd*G;
     (p.kdif)*a;
     p.kpd*Ra;
     p.kin;
     p0*p.kon*(1+p.f*(Ra/p.K)^2)/(1+(Ra/p.K)^2);
     p0ON*p.koff]; %*((1+(Ra/p.K)^2)/(1+p.f*(Ra/p.K)^2))];
end