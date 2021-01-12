%% Simulations of pendulum test using a muscle model 
%  Jente Willaert - 21/10/2020
%  clear all; close all; clc; 
path = 'C:\Users\u0125183\Documents\MATLAB\SimPendulum\'
addpath('C:\Users\u0125183\Documents\MATLAB\SimPendulum\BK/')
%% Input 
s.nu = 'TD5';   % subject number/ name 
s.tr = [1];     % subject trials (number of trials)
opt  = '';      % Option used as name to save results 
w1   = 0.7;
w2   = 0.7;     % Weigths in cost function

params = ImportParameters(s.nu);    % Input parameters (mtot,lc, l, age, m, RG, SE, Nmr, z)

for j = 1:length(s.tr)
    
    %% Prepare experimental data for simulation
    % Load experimental data + SRS on or off
    [q_exp_r, t_exp, t_span, on_srs] = LoadExpData(s.nu,s.tr,j,0, params, path);      % 1 if you want to plot experimental data
    
    % Discretised time: interpolate experimental data at discr. time using spline 
    [q_exp, qdot_exp,N, tvect, dt] = ExpDatAtDiscrTime(t_span,t_exp,q_exp_r);       % dt = 0.005
    
    % Define phases of pendulum (initial state, end of first swing) 
    [x0, N_1] = PendulumPhases(q_exp, qdot_exp, N);

    %% Formulate OCP 
    import casadi.*;        % Import casadi libraries
    opti = casadi.Opti();   % Initialise opti structure
    
    %% States using torques 
    
    % States - mesh points
    q     = opti.variable(1,N);
    qd    = opti.variable(1,N);
    T_srs = opti.variable(1,N); 
    Tq    = opti.variable(1,N); 
    Tdq   = opti.variable(1,N);
    Xe    = opti.variable(1,N);
    Xedot = opti.variable(1,N);
    
    % States - collocation points
    
    % Controls
    Tbq   = opti.variable(1); 
    Rq    = opti.variable(1);
    dRq   = opti.variable(1);
    Kr1   = opti.variable(1);
    Kr2   = opti.variable(1);
    K     = opti.variable(1); 
    C     = opti.variable(1);
    
    % Bounds
    opti.subject_to(-4*pi < q   < 4*pi);  
    opti.subject_to(-300  < qd  < 300);
    opti.subject_to(-4*pi < Xe  < 4*pi); 
    opti.subject_to(-300  < Xedot < 300);
    opti.subject_to(0.01  < Tbq < 10);
    opti.subject_to(0.01  < Rq  < 10);
    opti.subject_to(0.01  < dRq < 10); 
    opti.subject_to(0.01  < K   < 100);
    opti.subject_to(0.01  < C   < 100);
    opti.subject_to(-160  < Kr1 < -110);
    opti.subject_to(-30   < Kr2 < -30); 
    
    % Bounds on initial states
    opti.subject_to(q(1)     == x0(1));
    opti.subject_to(qd(1)    == x0(2));
    opti.subject_to(Xe(1)    == 0); 
    opti.subject_to(Xedot(1) == 0);
    
    % Initial guess
    opti.set_initial(q, q_exp);
    opti.set_initial(qd,qdot_exp);
    opti.set_initial(Tbq, 0);
    opti.set_initial(Rq, 0);
    opti.set_initial(dRq,0);
    opti.set_initial(K,0.01);
    opti.set_initial(C,1);
    opti.set_initial(Kr1, -160);
    opti.set_initial(Kr2, 10);
    
    %% Define problem 
    
    % Maxwell damping
    Xedot  = -(K/C)*Xe + qd;
    F_maxw  = K*Xe;
    x_v_dot = F_maxw/C;
    x_v     = q - Xe;
    opti.subject_to(Xedot(1:N-1)* dt + Xe(1:N-1) == Xe(2:N));
        
    % Tlim - limit torque to prevent hyperextension / hyperflexion of the leg
    knee_r_range = [Kr1 Kr2]*pi/180;
    klim         = 3;
    theta_ref    = (knee_r_range(2) - knee_r_range(1))/2;
        
    qrel = q - (knee_r_range(1) + knee_r_range(2))/2;
    Tlim = -exp(klim*(qrel-theta_ref)) + exp(klim*(-qrel-theta_ref));
        
    % Tsrs
    qc      = -1.5*pi/180;         % Qc = 1.5 graad = critical angle (Tsrs = ct*delta q if stretch is smaller than Qc / Tsrs = ct * Qc if stretch is larger than Qc)
    delta_q = q - q(1);
    kSRS    = 0.67*180/pi;         % ksrs = 0.67 nm/graad
              
        % Fase 1
        for k = 1:N_1
            opti.subject_to(T_srs(k)  == -((0.5*tanh(10*(delta_q(k)-qc))+0.5)*Tbq*kSRS*delta_q(k) + (0.5*tanh(10*(qc-delta_q(k)))+0.5)*kSRS*qc*Tbq));
            opti.subject_to(qd(k) < 0);
        end
        
        % Fase 2 - exponential decay of srs
        for k = N_1:N-1
            T_srs_dot = -T_srs(k)/0.05;
            opti.subject_to(T_srs(k+1) == T_srs(k) + T_srs_dot*dt);
        end
        
    % T reflex (only quadriceps)
    Trefl = [];
    TRH   = [];
        
        for k = 1: N-1
            opti.subject_to(Tq(k) == on_srs * T_srs(k) - F_maxw(k));
        end
        for k = 2:N                 % Toegevoegd op 14-4-2020
            opti.subject_to(Tq(k) == Tq(k-1) + Tdq(k)*dt);
        end
        
        for k = 1:N
            if k > 16
                Tq_delayed  = Tq(k-16);
                Tdq_delayed = Tdq(k-16);
            else
                Tq_delayed = 0;
                Tdq_delayed = 0;
            end
            
            Treflex = Rq* Tq_delayed + dRq*Tdq_delayed ;       % Reflex torque = force (srs + damping ) + derivative of force
            Trq  = (0.5*tanh(10*Treflex)+0.5) * Treflex;       % only positive torques
            Trefl= [Trefl; Trq];
        end
        
        % I adapted to optimized mass
        I = params.m*params.RG*params.RG + params.m*params.lc*params.lc;                   % Inertia of lower limb + foot
        
        % IP dynamics
        qdd = (-params.m*params.g*params.lc*cos(q))/I + Tbq/I - F_maxw/I + Tlim/I + Trefl'/I + on_srs * T_srs/I;
        
        % backward euler
        opti.subject_to(qd(1:N-1)*dt +q(1:N-1) == q(2:N));
        opti.subject_to(qdd(1:N-1)*dt +qd(1:N-1) == qd(2:N));
        
        % objective function
        qerror       = q - q_exp;
        qdoterror    = qd - qdot_exp;
        qf1_error    = q(N_1) - q_exp(N_1);                 % first swing
        qf2_error    = q(end-100:end) - q_exp(end-100:end); % resting angle
        qdotf1_error = qd(N_1);
        
        s1           = sumsqr(qdot_exp) / sumsqr(q_exp);    % scaling factor
        J1           = w1* sumsqr(qerror) + (1-w1) * sumsqr(qdoterror)/s1;
        J2           = (sumsqr(qf1_error)*length(q_exp) + sumsqr(qf2_error)*length(q_exp)/100);
        J            = w2* J1 + (1-w2) * J2 + 0.01*Tbq + 0.1*sumsqr(Trefl)/length(q_exp);
        opti.minimize(J);

        % options for IPOPT
        options.ipopt.tol = 1*10^(-10);
        options.ipopt.linear_solver = 'mumps';
        options.ipopt.linear_solver = 'ma57';
        % options.ipopt.hessian_approximation = 'limited-memory';
        
        % Solve the OCP
        opti.solver('ipopt',options);
        sol = opti.solve();
        
        %% Results
        
        qsim   = sol.value(q)';    qdsim  = sol.value(qd)';
        Xe     = sol.value(Xe);    Xe_dot = sol.value(Xedot); 
        tsim   = tvect;            Tbq    = sol.value(Tbq);
        Rq     = sol.value(Rq);    dRq    = sol.value(dRq);
        K      = sol.value(K);     C      = sol.value(C);
        Kr1    = sol.value(Kr1);   Kr2    = sol.value(Kr2);
        F_maxw = sol.value(F_maxw);
        
        % Torques based on optimized parameters
        Trefl  = sol.value(Trefl);  Tsrs   = sol.value(T_srs);
        Tlim   = sol.value(Tlim);
                
        figure()
        plot(tsim,q_exp*180/pi, 'k','LineWidth',1.5); hold on;
        plot(tsim,qsim*180/pi,'r','LineWidth',1.5); hold on;
        legend({'Q: Exp','Q: Sim'})
        ylabel('Angle [{\circ}]');   
end
