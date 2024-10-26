%This code determined the QPD sensitivity and trap stiffness of a trapped
%microsphere in water
%This code was written in Rheology and Light scattering lab in Raman
%Research Institute, India on 26-10-2024.
%Procedure adapted from 'Optical trapping' by Keir C. Neuman and Steven M.
%Block (https://pubs.aip.org/aip/rsi/article/75/9/2787/351902/Optical-trapping)
clear;
close all;

%input parameters
R = (1/2)*1e-6; %Probe radius [m]
T = 24; %bath temperature [*C]
T = T + 273.16; % K
kB = 1.38064852e-23; %Boltzman constant
eta = 0.89*1e-3; % viscosity of water [Pa-s]
samplingrate = 30000; % Data acquisition rate [Hz]
dat = load('data.dat'); %[Vx, Vy, Vsum] fluctuations of trapped microsphere


gamma=6*pi*eta*R; % viscous damping factor
D0=kB*T/(gamma); 
dt = 1/samplingrate; 
% Determining QPD sensitivity and trap stiffness along the x-direction
Vx = dat(:,1);
Vsum = dat(:,3);
x = Vx./Vsum;
x = x - repmat(mean(x),size(x,1),1); %Centering data to zero.
N = length(x);
t1 = (0:1/samplingrate:(N/samplingrate))';

% Calculating the power spectral density
ll = 8;
NN = N/(ll*samplingrate);
L = ll*samplingrate;
xpsd = zeros(L,NN);
for i = 1:NN
    xx(:,i) = x(1+(i-1)*L:i*L);
    X = fft(xx(:,i));
    % xpsd(:,i) = ((X.*conj(X)).*dt.^2)./(N*dt);
    [Pxx(:,i),ff1] = periodogram(xx(:,i),[],length(xx(:,i)),samplingrate);
end
Pxx1 = mean(Pxx')';
% xpsd1 = mean(xpsd')';
loglog(ff1,Pxx1,'.-');
xlabel('f [Hz]','FontSize',16);
ylabel('PSD [QPD unit^2/Hz]','FontSize',16);
id1 = find(ff1 == 10);
id2 = find(ff1 == 11000);
ff2 = ff1(id1:id2+1,1);
Pxx2 = Pxx1(id1:id2+1,1);
loglog(ff2,Pxx2,'.-');
xlabel('f [Hz]','FontSize',16);
ylabel('PSD [AU^2/Hz]','FontSize',16);
%Determining QPD caliberation factor
[Pv2,cal_factor2] = Power(ff2,Pxx2,gamma,kB,T);
%Determining trap stiffness using PSD method
[D2,k_psd] = lorentz(ff2,Pxx2,gamma);
function [D,kappa] = lorentz(f3,xpsd4,gamma)
plot(f3,xpsd4,'.-','Linewidth',0.5);
set(gca,'XScale','Log');
set(gca,'YScale','Log');
hold on;
xpsd4 = xpsd4.*1e8;
%Lorentzian fitting
a = 100;
b = 1500;
clear('xData','yData');
[xData, yData] = prepareCurveData(f3, xpsd4);
ft = fittype( 'a/(b^2 + x^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = [a b];
[fitresult, gof] = fit(xData, yData, ft, opts);
fc2 = fitresult.b;
D2 = fitresult.a;
D2 = D2/1e8;
% d = fitresult.d;
% d = d/1e8;
out = D2./(fc2^2 + f3.^2); %+ d;
plot(f3,out,'k','LineWidth',3);
xlabel('f [Hz]');
ylabel('S_{xx} [nm^2/Hz]');
hold off;
xpsd4 = xpsd4./1e8;
D = pi^2*D2;
kappa = 2*pi*gamma*fc2; %N/m
XX = ['k =',num2str(kappa*1e6),' pN/um.'];
disp('...')
disp('PSD analysis')
disp(XX);
end

function [Pv,cal_factor] = Power(f1,xpsd2,gamma,kB,T)
f2xpsd = (f1.^2).*(xpsd2);
% inum = 20000;
% fnum = 35000;
% f_plateu = f1(inum:fnum);
% y_plateu = f2xpsd(inum:fnum);
loglog(f1,f2xpsd,'.-');
x = [20000; 35000];
%f_plateu = f1(ceil(x(1,1)):floor(x(2,1)));
y_plateu = f2xpsd(ceil(x(1,1)):floor(x(2,1)));
Pv = mean(y_plateu);
%QPD Caliberation factor
beta = sqrt((Pv*pi^2*gamma)/(kB*T)); % Volt/m
%beta = beta; 
cal_factor = (1/beta);
XX = ['cal factor = ',num2str(cal_factor),' m/V.'];
disp(XX);
end

