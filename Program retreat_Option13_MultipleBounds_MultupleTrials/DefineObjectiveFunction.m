function J = DefineObjectiveFunction(x,xd,qspline, qdspline, info, vMtilda, a_ext, a_flex, kR)
%Define objective function that should be minimized 
w_q  = info.wq;
w_qd = info.wqd;

q_error  = x - qspline'; 
qd_error = xd - qdspline'; 

J = w_q * sumsqr(q_error) + w_qd * sumsqr(qd_error) + 0.001 * (sumsqr(vMtilda)) + 0.1*(a_ext + a_flex + kR); 



end

