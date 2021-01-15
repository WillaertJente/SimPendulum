a = 0.2;
FMo = 1000;
FMltilda = 1;
kSRS = 280;

deltaL = -0.1:0.001:0.7;
Fsrs = zeros(size(deltaL));

b = 1000;

for i = 1:length(deltaL)
dLm = deltaL(i);
Fsrs(i)     =(0.5*tanh(b*(-dLm+5.7*10^(-3)))+0.5)*dLm.*FMltilda*a*kSRS + (0.5*tanh(b*(dLm - 5.7*10^(-3)))+0.5)*5.7*10^(-3)*a*FMltilda*kSRS;

end

figure()
plot(deltaL, Fsrs)
xlabel('\Delta L [-]')
ylabel('F_{SRS}')