function [x0, N_1] = PendulumPhases(q_exp, qdot_exp, N, bool)
%Define initial state of pendulum and end state of phase 1
%   1. Define initial state
%   2. Define end of phase 1

% Initial state of pendulum
x0 = [q_exp(1) qdot_exp(1)]; 

% Define end of first phase (for SRS)
for i = 1:N
    if qdot_exp(i)*qdot_exp(i+1) < 0 && i > 20
        break;
    end
end

N_1 = i

% Plot phase
if bool == 1
    figure()
    plot(q_exp,'k','LineWidth',1.5)
    hold on
    line([i i],[-2.5 0],'Color','r','LineWidth',1.5)

end

