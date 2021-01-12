clear all; close all;

Del1 = importdata('Del1_dMLMT')

shoulder_elev = Del1.data(:,2);
dM_Del1 = Del1.data(:,3);
LMT_Del1 = Del1.data(:,4);

figure(3)
subplot(2,1,1)
plot(shoulder_elev, dM_Del1); hold on;
ylim([-0.1 0.1])

title('Moment arm')
subplot(2,1,2)
plot(shoulder_elev, LMT_Del1); hold on;
title('Muscle tendon length')

x_Del1= muscleFitCosine(shoulder_elev, dM_Del1, LMT_Del1);


theta = pi/180*shoulder_elev;

[dM_est,LMT_est] = evaluate_dM_LMT(x_Del1,theta);

subplot(2,1,1)
plot(shoulder_elev, dM_est); hold on;
ylim([-0.1 0.1])

title('Moment arm')
subplot(2,1,2)
plot(shoulder_elev, LMT_est); hold on;
title('Muscle tendon length')



Del3 = importdata('Del3_dMLMT')

shoulder_elev = Del3.data(:,2);
dM_Del3 = Del3.data(:,4);
LMT_Del3 = Del3.data(:,3);

figure(4)
subplot(2,1,1)
plot(shoulder_elev, dM_Del3); hold on;
ylim([-0.1 0.1])

title('Moment arm')
subplot(2,1,2)
plot(shoulder_elev, LMT_Del3); hold on;
title('Muscle tendon length')

x_Del3= muscleFitCosine(shoulder_elev, dM_Del3, LMT_Del3);

theta = pi/180*shoulder_elev;

[dM_est,LMT_est] = evaluate_dM_LMT(x_Del3,theta);

subplot(2,1,1)
plot(shoulder_elev, dM_est); hold on;
title('Moment arm')
ylim([-0.1 0.1])
subplot(2,1,2)
plot(shoulder_elev, LMT_est); hold on;
title('Muscle tendon length')


TriLat = importdata('TriLat_dMLMT')

elbow_flexion = TriLat.data(:,2);
dM_TriLat = TriLat.data(:,3);
LMT_TriLat = TriLat.data(:,4);

figure(1)
subplot(2,1,1)
plot(elbow_flexion, dM_TriLat); hold on;
ylim([-0.1 0.1])

title('Moment arm')
subplot(2,1,2)
plot(elbow_flexion, LMT_TriLat); hold on;
title('Muscle tendon length')

x= muscleFitCosine(elbow_flexion, dM_TriLat, LMT_TriLat);


theta = pi/180*elbow_flexion;
[dM_est,LMT_est] = evaluate_dM_LMT(x,theta);


subplot(2,1,1)
plot(elbow_flexion, dM_est); hold on;
ylim([-0.1 0.1])

title('Moment arm')
subplot(2,1,2)
plot(elbow_flexion, LMT_est); hold on;
title('Muscle tendon length')



BicLong = importdata('BicLong_dMLMT')

elbow_flexion = BicLong.data(:,2);
dM_BicLong = BicLong.data(:,4);
LMT_BicLong = BicLong.data(:,3);

figure(2)
subplot(2,1,1)
plot(elbow_flexion, dM_BicLong); hold on;
ylim([-0.1 0.1])

title('Moment arm')
subplot(2,1,2)
plot(elbow_flexion, LMT_BicLong); hold on;
title('Muscle tendon length')

x_BicLong= muscleFitCosine(elbow_flexion, dM_BicLong, LMT_BicLong);


theta = pi/180*elbow_flexion;
[dM_est,LMT_est] = evaluate_dM_LMT(x_BicLong,theta);


subplot(2,1,1)
plot(elbow_flexion, dM_est); hold on;
title('Moment arm')
subplot(2,1,2)
plot(elbow_flexion, LMT_est); hold on;
title('Muscle tendon length')






dM_LMT_params = [ x_Del1'; x_Del3';x'; x_BicLong';];
save('dM_LMT_params.mat','dM_LMT_params');




figure(5)
VMT = evaluate_VMT(x_Del3',theta,1);
plot(theta,VMT);
