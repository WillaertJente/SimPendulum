function cost = getFitCost(x,theta,dM,LMT)


[dM_est,LMT_est] = evaluate_dM_LMT(x,theta);
cost = sumsqr(2*(dM - dM_est)) + sumsqr((LMT - LMT_est));

end