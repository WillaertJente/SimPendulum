function dadt = ActivationOde(e,a)

b = 10;

tact = 0.015;
tdeact = 0.06;

d1 = 1./(tact*(0.5+1.5*a));
d2 = (0.5+1.5*a)/tdeact;
f = 0.5*tanh(b*(e-a));
dadt = (d1.*(f+0.5) + d2.*(-f+0.5)).*(e-a);