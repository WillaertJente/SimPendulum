function [q,qd,lMtilda] = forwardSim(q_init,qd_init,lMtilda_init,kFpe,a,data_exp, coeff_LMT_ma, params_OS, shift, dt, N, B)
nMuscles = 2;
nJoints = 1;
nStates = nMuscles + 2*nJoints;
a_ext = a(1); 
a_flex = a(2); 
import casadi.*;


q = NaN(nJoints,N+1); q(:,1) = q_init;
qd = NaN(nJoints,N+1); qd(:,1) = qd_init;
lMtilda = NaN(nMuscles,N+1); lMtilda(:,1) = lMtilda_init;
guess = 0.1*ones(3*nStates+2*nMuscles,1);
for i = 1:N
    
    Urf = MX.sym('Urf',3*nStates+2*nMuscles);

    % Unknown state derivatives at mesh interval start
    qd_implicit = Urf(1:nJoints);
    qdd_implicit = Urf(nJoints+1:2*nJoints);
    vMtilda_implicit = Urf(2*nJoints+1:2*nJoints+nMuscles);

    % Unknown states at mesh interval end
    q_next = Urf(2*nJoints+nMuscles+1:3*nJoints+nMuscles);
    qd_next = Urf(2*nJoints+nMuscles+nJoints+1:4*nJoints+nMuscles);
    lMtilda_next = Urf(2*nJoints+nMuscles+2*nJoints+1:4*nJoints+2*nMuscles);

    % Unknown state derivatives at mesh interval end
    qd_implicit_next = Urf(4*nJoints+2*nMuscles+1:4*nJoints+2*nMuscles+nJoints);
    qdd_implicit_next = Urf(4*nJoints+2*nMuscles+nJoints+1:4*nJoints+2*nMuscles+2*nJoints);
    vMtilda_implicit_next = Urf(4*nJoints+2*nMuscles+2*nJoints+1:4*nJoints+2*nMuscles+2*nJoints+nMuscles);

    lM_projected = Urf(4*nJoints+2*nMuscles+2*nJoints+nMuscles+1:4*nJoints+2*nMuscles+2*nJoints+2*nMuscles);
    lM_projected_next = Urf(4*nJoints+2*nMuscles+2*nJoints+2*nMuscles+1:4*nJoints+2*nMuscles+2*nJoints+3*nMuscles);

    
    errorDyn_meshStart = CalculateMusculoSkeletalDynamics(q(:,i),qd_implicit,qdd_implicit,lMtilda(:,i),lM_projected,kFpe, vMtilda_implicit,a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B);
    errorDyn_meshEnd = CalculateMusculoSkeletalDynamics(q_next,qd_next,qdd_implicit_next,lMtilda_next,lM_projected_next,kFpe, vMtilda_implicit_next,a_ext, a_flex, data_exp, coeff_LMT_ma, params_OS, shift, B);
    errorVel_meshStart = qd_implicit - qd(:,i);
    errorVel_meshEnd = qd_implicit_next - qd_next;
    dlMdt_implicit = CalculateDLMDT(vMtilda_implicit, params_OS); 
    dlMdt_implicit_next = CalculateDLMDT(vMtilda_implicit_next, params_OS); 
    errorIntegration = [(qd_implicit + qd_implicit_next)*dt/2 + q(:,i) - q_next;
                        (qdd_implicit + qdd_implicit_next)*dt/2 + qd(:,i) - qd_next;
                        (dlMdt_implicit + dlMdt_implicit_next)*dt/2 + lMtilda(:,i) - lMtilda_next];
                        
    lMo    = params_OS.MT(2,:); 
    
    rf = rootfinder('rf','kinsol',struct('x',Urf,'g',[errorDyn_meshStart'; errorDyn_meshEnd'; errorIntegration;errorVel_meshStart;errorVel_meshEnd]),struct('abstol',1e-14,'error_on_fail',false));
    
        
    guess = [qd(:,i);       0;            zeros(nMuscles,1); 
             q(:,i) ; qd(:,i); lMtilda(:,i);
             qd(:,i); 0; zeros(nMuscles,1);
             lMtilda(:,i).*lMo(1);  lMtilda(:,i).*lMo(1)];

    solution = full(rf(guess,[]));
    q(:,i+1) = solution(2*nJoints+nMuscles+1:3*nJoints+nMuscles);
    qd(:,i+1) = solution(2*nJoints+nMuscles+nJoints+1:4*nJoints+nMuscles);
    lMtilda(:,i+1) = solution(2*nJoints+nMuscles+2*nJoints+1:4*nJoints+2*nMuscles);
    guess = 0.1*ones(3*nStates+2*nMuscles,1);
end
end