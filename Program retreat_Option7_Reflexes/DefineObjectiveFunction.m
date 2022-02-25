function J = DefineObjectiveFunction(x,xd,data_exp, info, vMtilda) 
%Define objective function that should be minimized 
w_q  = info.wq;
w_qd = info.wqd;

q_error  = x - data_exp.qspline'; 
qd_error = xd - data_exp.qdspline'; 

J = w_q * sumsqr(q_error) + w_qd * sumsqr(qd_error) + 0.001 * (sumsqr(vMtilda)) ; 



end

