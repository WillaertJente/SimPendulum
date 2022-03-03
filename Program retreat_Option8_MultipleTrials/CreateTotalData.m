function [data_exp] = CreateTotalData(data_p1, data_p2)
%UNTITLED Summary of this function goes here
data_exp.qspline  = [data_p1.qspline; data_p2.qspline];
data_exp.qdspline = [data_p1.qdspline; data_p2.qdspline];
data_exp.Nspline  = [data_p1.Nspline data_p2.Nspline];
data_exp.offset   = [data_p1.offset; data_p2.offset];
data_exp.x0       = [data_p1.x0; data_p2.x0]; 
data_exp.N_1      = [data_p1.N_1 data_p2.N_1];
end

