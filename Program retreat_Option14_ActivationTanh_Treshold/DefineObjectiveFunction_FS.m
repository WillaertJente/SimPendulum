function J = DefineObjectiveFunction_FS(x,xd,data_exp, info, vMtilda, a_ext, a_flex, kR)
%Define objective function that should be minimized 
w_q  = info.wq;
w_qd = info.wqd;

q_error  = x - data_exp.qspline'; 
qd_error = xd - data_exp.qdspline'; 
FS_error = data_exp.qspline(data_exp.N_1) - x(data_exp.N_1); 

J = w_q * sumsqr(q_error) + w_qd * sumsqr(qd_error) + 0.001 * (sumsqr(vMtilda)) + 0.1*(a_ext + a_flex + kR) + 10*FS_error; 



end

