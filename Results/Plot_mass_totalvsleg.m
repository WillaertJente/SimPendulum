%% Plot mass total vs mass leg 
input = [0.1735 17.60; 0.4143 25.2; 0.7196 30.10; 1.0760 33.80; 1.6756 34; 1.3250 38; 1.8066 52.4; ];

figure()
subplot(221)
scatter(input(:,1),input(:,2),'Filled') % outlier = TD
ylabel('total mass'); xlabel('leg mass')
%% Plot mass total vs damping
damping = [0.0532; 0.0730; 0.0112; 0.0307; 0.0126; 0.1626;0.0234];

subplot(222)
scatter(input(:,2),damping,'Filled')
xlabel('total mass'); ylabel('damping')

%% plot mass leg vs. damping
subplot(223)
scatter(input(:,1),damping,'Filled')
xlabel('Leg mass'); ylabel('damping')

%% Ratios
% input = mass total leg sorted from small to large 
figure()
subplot(221)
factor = damping./input(:,2)
scatter(1:length(input), factor,'Filled')
title('Total mass/ damping')

subplot(222)
input = [0.1735 17.60; 0.4143 25.2; 0.7196 30.10; 1.0760 33.80;1.3250 38; 1.6756 34; 1.8066 52.4; ];
factor = damping./input(:,1)
scatter(1:length(input), factor,'Filled')
title('Leg mass / damping')

%% Plot leg mass tov leg inertia 
input = [1.6756 0.0758; 1.0760 0.0389; 0.1735 0.0029; 1.3250 0.0509; 0.4143 0.0096; 1.8066 0.0689; 0.7196 0.0215 ]
figure(10)
subplot(221)
scatter(input(:,1),input(:,2),'Filled'); hold on
ylabel('Inertia'); xlabel('Leg mass')

subplot(222)
scatter(input(:,2),damping,'Filled'); hold on
xlabel('Inertia'); ylabel('Damping') 