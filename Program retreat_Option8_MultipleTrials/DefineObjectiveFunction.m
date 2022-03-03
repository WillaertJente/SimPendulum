function J = DefineObjectiveFunction(x,xd,qspline, qdspline, info, vMtilda) 
%Define objective function that should be minimized 
w_q  = info.wq;
w_qd = info.wqd;

q_error  = x - qspline'; 
qd_error = xd - qdspline'; 

J = w_q * sumsqr(q_error) + w_qd * sumsqr(qd_error) + 0.001 * (sumsqr(vMtilda)) ; 



end

