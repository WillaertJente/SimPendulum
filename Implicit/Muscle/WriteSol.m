function [Results] = WriteSol(result,N_all,N_1_all,tr)
%Function to write results 

for j = 1:length(tr)
    if j == 1
        start(j)     = 1;
        stop(j)      = N_all(j);
        N_1_start(j) = 1;
        N_1_stop(j)  = N_all(j)-N_1_all(j);
    else
        start(j)     = stop(j-1)+1;
        stop(j)      = start(j)-1 + N_all(j);
        N_1_start(j) = N_1_stop(j-1)+1;
        N_1_stop(j)  = N_1_start(j)-1 + N_all(j)-N_1_all(j);
    end

%% Results trial 1
if j == 1 
    Results.x1                = result(start(j):stop(j));
    Results.xd1               = result(stop(j)+1:2*stop(j)); 
    Results.lMtilda_ext1      = result(2*stop(j)+1:3*stop(j));
    Results.lMtilda_flex1     = result(3*stop(j)+1:4*stop(j));
    Results.Fsrs21            = result(4*stop(j)+1:4*stop(j)+stop(j)-N_1_stop(j));
    Results.Fsrs_d1           = result(5*stop(j)-N_1_stop(j)+1:6*stop(j)-N_1_stop(j));
    Results.a_ext1            = result(6*stop(j)-N_1_stop(j)+1:7*stop(j)-N_1_stop(j));
    Results.lM_projected_ext1 = result(7*stop(j)-N_1_stop(j)+1:8*stop(j)-N_1_stop(j));
    Results.lM_projected_flex1= result(8*stop(j)-N_1_stop(j)+1:9*stop(j)-N_1_stop(j));
    Results.act1              = result(9*stop(j)-N_1_stop(j)+1:10*stop(j)-N_1_stop(j));
    Results.dt11              = result(10*stop(j)-N_1_stop(j)+1);
    Results.vMtilda_ext1      = result(10*stop(j)-N_1_stop(j)+2:11*stop(j)-N_1_stop(j)+1);
    Results.vMtilda_flex1     = result(11*stop(j)-N_1_stop(j)+2:12*stop(j)-N_1_stop(j)+1);
    Results.a_ext_01          = result(12*stop(j)-N_1_stop(j)+2);
    Results.a_flex1           = result(12*stop(j)-N_1_stop(j)+3);
    Results.kFpe_ext1         = result(12*stop(j)-N_1_stop(j)+4);
    Results.kFpe_flex1        = result(12*stop(j)-N_1_stop(j)+5); 
    Results.Rk1               = result(12*stop(j)-N_1_stop(j)+6);
    





end
end


