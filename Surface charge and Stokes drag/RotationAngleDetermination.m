clear;
close all;

R = (1/2)*1e-6; %Probe radius [m]
fAC = 7991.9; %AC signal frequency
T = 24; %bath temperature
T = T + 273.16;
kB = 1.38064852e-23; %Boltzman constant
samplingrate = 30000;
dat = load('trj.dat');
%% 

dt = 1/samplingrate;
Vx = dat(:,2); %<--y [wrt AC field]
Vy = dat(:,1); %<--x [wrt AC field]
Vsum = dat(:,3);
x = (Vx./Vsum); % In QPD units
y = (Vy./Vsum); % In QPD units
%centering data to zero
x = x - repmat(mean(x),size(x,1),1);
y = y - repmat(mean(y),size(y,1),1);
data_centered = [x y];
%
N = length(x);
ll = 10; % Choose it such that it is divisible by dt*N = 80
NN = floor(N/(ll*samplingrate));
L = ll*samplingrate;
%% Determine rotation angle
rotationangle = 0; %Can be any initial value
rotation_matrix = [cos(-rotationangle), -sin(-rotationangle); sin(-rotationangle), cos(-rotationangle)];
rotated_data = data_centered * rotation_matrix;
theta = -0.005;
rotation_matrix = [cos(-theta), -sin(-theta); sin(-theta), cos(-theta)];
temp = 0;
kk = 1;
ll = 1;
while ll == 1
    tempdata = rotated_data(:,1);
    tempdata2 = rotated_data(:,2);
    for i = 1:NN
        xx1 = tempdata(1+(i-1)*L:i*L);
        yy1 = tempdata2(1+(i-1)*L:i*L);
        [yPxx(:,i),yff1] = periodogram(xx1,[],length(xx1),samplingrate);
        [xPyy(:,i),yff1] = periodogram(yy1,[],length(yy1),samplingrate);
    end
    ypsd1 = mean(yPxx')';
    xpsd1 = mean(xPyy')';
    id1_1 = find(yff1 == 10);
    id2_1 = find(yff1 == 9600);
    r_ypsd1 = ypsd1(id1_1:id2_1,1);
    r_xpsd1 = xpsd1(id1_1:id2_1,1);
    r_yff1 = yff1(id1_1:id2_1,1);
    [~, fid] = min(abs(r_yff1 - fAC));

    figure(1);
    subplot(1,3,1);
    scatter(rotated_data(1:10000,1),rotated_data(1:10000,2));
    xlabel('x [QPD units]');
    ylabel('y [QPD units]');
    title('Position fluctuation');
    axis equal;
    box on;
    subplot(1,3,2);
    loglog(r_yff1,r_ypsd1,'o-');
    xlabel('f [Hz]');
    ylabel('S_{yy} [QPD units^2/Hz]');
    title('y-direction');
    box on;
    subplot(1,3,3);
    loglog(r_yff1,r_xpsd1,'o-');
    xlabel('f [Hz]');
    ylabel('S_{xx} [QPD units^2/Hz]');
    title('x-direction');
    box on;
    shg;
    temp2 = temp;
    temp = r_ypsd1(fid,1); %sum(r_ypsd1(fid:fid,1)); 
    fprintf('\n temp = %f, Rotation angle = %f \t \n', abs(temp*1e9),abs(rotationangle));
    if kk>2
        if temp2*1e9 < temp*1e9
            ll = 0;
            break;
        end
    end
     rotated_data = rotated_data * rotation_matrix;
     rotationangle = rotationangle + theta;
     kk = kk+1;
     clear('yPxx','r_ypsd1');
end
fprintf('Final rotation angle = %f \n', abs(rotationangle - theta));
