%% SCript to save BK data for simulations
map   = 'C:\Users\u0125183\Box\PhD 1\Simulations Pendulum Test\CP8/'
nr_tr = [2 3 4 5];  

for i = nr_tr
%% Import BK data 
BK_name = [map, 'BK_Trial0',num2str(i),'_BodyKinematics_pos_global.sto']
BK      = importdata(BK_name);

%% Plot BK data
col_tibia = find(strcmp('tibia_r_Oz',BK.colheaders))

figure(i)
plot(BK.data(:,col_tibia))

data = [BK.data(:,1) BK.data(:,col_tibia)];

end
