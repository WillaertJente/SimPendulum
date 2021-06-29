
Fpe = 0:0.01:2;
Fpe_lim = (- 0.5*tanh(10*(Fpe-1.5)) + 0.5).*Fpe + (0.5*tanh(10*(Fpe-1.5)) + 0.5)*1.5;

figure()
plot(Fpe, 'b'); hold on;
plot(Fpe_lim, 'k')

a_ext_0 = 0.228;
dLm = [0: 0.0001: 0.05];
FMltilda_ext = 1;
kSRS = 280;

Fsrs     =(0.5*tanh(1000*(-dLm+5.7*10^(-3)))+0.5).*dLm.*FMltilda_ext*a_ext_0*kSRS + ...
    (0.5*tanh(1000*(dLm - 5.7*10^(-3)))+0.5)*5.7*10^(-3).*a_ext_0.*FMltilda_ext*kSRS;

figure()
plot(Fsrs)