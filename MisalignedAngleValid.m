% Validate misalignment angle and lever arm of reference IMU to target IMU

%%
clear all;
close all;
clc;

%% file load and data interpolation
[A15_file, A15_path] = uigetfile('*.txt');
A15_file = [A15_path, A15_file];
[tar_file, tar_path] = uigetfile('*.txt');
tar_file = [tar_path, tar_file];

A15 = importdata(A15_file);
[~,uni_id] = unique(A15(:,1));
A15 = A15(uni_id,:);
A15(:,10) = smooth_angle(A15(:,10));

target = importdata(tar_file);
[~,uni_id] = unique(target(:,1));
target = target(uni_id,:);

id_start = 0;
id_end = 0;
if target(1,1) < A15(1,1)
    id_start = find(target(:,1) > A15(1,1));
    target = target(id_start(1):end,:);
end

if target(end,1) > A15(end,1)
    id_end = find(target(:,1) < A15(end,1));
    target = target(1:id_end(end),:);
end
target(:,10) = smooth_angle(target(:,10));

tmp = zeros(size(target));
tmp(:,1) = target(:,1);
for i = 2:size(target,2)
    tmp(:,i) = interp1(A15(:,1),A15(:,i),target(:,1));
end
A15 = tmp;

%% convert to NED coordinate
a = 6378137.0;
b = 6356752.3141;
e2 = 0.00669437999013;

phi1 = A15(:, 2) * pi / 180.0;
theta1 = A15(:, 3) * pi / 180.0;
h1 = A15(:, 4);
Rn1 = a ./ sqrt(1 - e2 * sin(phi1) .* sin(phi1));
A15(:, 2) = (Rn1 + h1) .* cos(phi1) .* cos(theta1);
A15(:, 3) = (Rn1 + h1) .* cos(phi1) .* sin(theta1);

phi2 = target(:, 2) * pi / 180.0;
theta2 = target(:, 3) * pi / 180.0;
h2 = target(:, 4);
Rn2 = a ./ sqrt(1 - e2 * sin(phi2) .* sin(phi2));
target(:, 2) = (Rn2 + h2) .* cos(phi2) .* cos(theta2);
target(:, 3) = (Rn2 + h2) .* cos(phi2) .* sin(theta2);

error = A15-target;
error(:,end) = wrap_angle(error(:,end)); % yaw error process
error = error(:,2:end);
mean_err = mean(error, 1);
disp('Mean position error (m):');
fprintf(sprintf('%0.3f, %0.3f, %0.3f\n', mean_err(1:3)));
disp('Mean velocity error (m/s):');
fprintf(sprintf('%0.3f, %0.3f, %0.3f\n', mean_err(4:6)));
disp('Mean attitude error (deg):');
fprintf(sprintf('%0.3f, %0.3f, %0.3f\n\n', mean_err(7:9)));
rms_err = rms(error,1);
disp('STD position error (m):');
fprintf(sprintf('%0.3f, %0.3f, %0.3f\n', rms_err(1:3)));
disp('STD velocity error (m/s):');
fprintf(sprintf('%0.3f, %0.3f, %0.3f\n', rms_err(4:6)));
disp('STD attitude error (deg):');
fprintf(sprintf('%0.3f, %0.3f, %0.3f\n', rms_err(7:9)));

figure;
plot(target(:,1),error(:,1:3));
grid on;
title('Pos error/m');
figure;
plot(target(:,1),error(:,4:6));
grid on;
title('Velocity error/(m/s)');
figure;
plot(target(:,1),error(:,7:9));
grid on;
legend('R','P','Y');
title('Attitude error/deg');
