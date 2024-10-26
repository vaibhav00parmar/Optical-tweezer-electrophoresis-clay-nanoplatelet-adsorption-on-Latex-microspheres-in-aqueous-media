%Written by vaibhav 16-09-2024
clear;
close all;
%
input_folder = 'D:\OT2024\Organised\Surface Charge Experiments\32 - 09092024\1 um PS in water\In field\Variable E\splitfiles5';
cd(input_folder);
%Experimental Parameters
rotationangle = 0.155;
d = 5.0; %Gap between the electrodes [mm]

k = 50.61*1e-6; % N/m 
R = (1/2)*1e-6; %Probe radius [m]
cal_factorx = 2.01615E-6; % [m/V] 
cal_factory = 2.39308E-6; % [m/V] 

V = 10.0;
E = (V/d)*1e3; %Applied Field [V/m]
fAC = 7991.9; %AC signal frequency %7991.9
samplingrate = 30000; %SamplingRate
T = 24; %bath temperature
T = T + 273.16;
kB = 1.38064852e-23; %Boltzman constant
%%
for kk = 1:1
    cd(input_folder);
    dat = load(sprintf('test_%d.lvm',kk));
    dt = 1/samplingrate;
    Vx = dat(:,2); %Parallel to the field %<--x [wrt AC field]
    Vy = dat(:,1); %Perpendicular to the appliued field %<--y [wrt AC field]
    Vsum = dat(:,3);
    x = (Vx./Vsum).*cal_factorx;
    y = (Vy./Vsum).*cal_factory;
    x = x - repmat(mean(x),size(x,1),1);
    y = y - repmat(mean(y),size(y,1),1);
    data_centered = [x y];
    %Rotate
    rotation_matrix = [cos(-rotationangle), -sin(-rotationangle); sin(-rotationangle), cos(-rotationangle)];
    rotated_data = data_centered * rotation_matrix;
    x = rotated_data(:,1);
    N = length(x);
    t1 = (0:1/samplingrate:(N/samplingrate))';
    t = t1(1:end-1,1);
    %%
    ll = 10; % Choose it such that it is divisible by dt*N = 80
    NN = N/(ll*samplingrate);
    L = ll*samplingrate;
    % xpsd = zeros(L,NN);
    for i = 1:NN
        xx(:,i) = x(1+(i-1)*L:i*L);
        [Pxx(:,i),ff1] = periodogram(xx(:,i),[],length(xx(:,i)),samplingrate);
    end
    %% mean power spectral density without moving average
    Pxx1 = mean(Pxx')';
    id1 = find(ff1 == 10);
    id2 = find(ff1 == 11000); %
    ff2 = ff1(id1:id2+1,1);
    Pxx2 = Pxx1(id1:id2+1,1);
    figure();
    loglog(ff2,Pxx2,'.-k');
    ylabel('S_{xx} [m^2/Hz]');
    xlabel('f [Hz]');
    hold on;
    %% If noise is present at certain frequencies [Else don't use this section]
    fid1=195;
    fid2=205;
    [d1,ix1] = min(abs(ff2-fid1));
    [d2,ix2] = min(abs(ff2-fid2));
    Pxx2(ix1:ix2,1) = NaN;
    ff2(ix1:ix2,1) = NaN;
    ff2 = rmmissing(ff2);
    Pxx2 = rmmissing(Pxx2);
    loglog(ff2,Pxx2,'o-b');
    %% Removing peak corresponding to f = fAC
    [d1,ix1] = min(abs(ff2-fAC));
    Pxx2_1 = Pxx2;
    for i = ix1-60:ix1+60
        Pxx2_1(i,1) = Pxx2(i-200,1);
    end
    loglog(ff2,Pxx2_1,'.b');
    %% Fitting lorentzian to determine fc and gamma (drag)
    Pxx2_2 = Pxx2_1.*1e18;
    custom_fit = @(a, b, x) a./(b^2 + x.^2);
    ft = fittype(custom_fit, 'independent', 'x', 'dependent', 'y');
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Lower = [0 0];
    opts.StartPoint = [10, 1100];
    [fitresult2, gof2] = fit(ff2, Pxx2_2, ft, opts);
    fc = fitresult2.b;
    D = fitresult2.a/1e18;
    %d = fitresult2.d/1e6;
    out = (D)./(fc^2 + ff2.^2);
    plot(ff2,out,'-b','LineWidth',2);
    hold off;

    gamma1 = (kB*T)/(pi^2*D); %<--Using this as drag coefficient
    gamma2 = k/(2*pi*fc);

    eta_eff1 = gamma1/(6*pi*R);
    eta_eff2 = gamma2/(6*pi*R);
    %% Calculations
    area1_1 = trapz(ff2,Pxx2_1);
    area2_1 = trapz(ff2,Pxx2);
    Pac = area2_1 - area1_1;
    Gammasq = (Pac*k)/(kB*T);
    Gamma = sqrt(Gammasq);
    Qeff = (Gamma/E)*(sqrt(2*kB*T*k*(1+(fAC/fc)^2)));
    e = 1.602176634E-19; %C
    Zeff = Qeff/e;

    res(kk,:) = [E fAC k Gammasq Zeff fc gamma1 gamma2 eta_eff1 eta_eff2];
    filename = sprintf('image_%d.png', kk);
    saveas(gcf, filename);
 end
