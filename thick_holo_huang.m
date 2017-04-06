% clear all;
% close all;


% This script is intended to calculate diffraction efficiency of thick hologram gratings  
%   from Huang and Gilbert's "Diffraction properties of substrate guided-wave holograms", 
% Opt. Eng, V.34, No10, 2891-2899, 1995

L = 647; % reference wavelength in nm
t = 7; % thickness in micron
n = 1.63; % avg index of refraction of recording medium
dn = 0.05; % index modulation
ar = 120; % value for Bragg angle in DEGREES % ??????????
ao = 0;

% Angular selectivity:
damin = 0; % angular deviation in degrees
damax = 5;
% Wavelength selectivity:
dLmin = 0; % wavelength deviation in nm
dLmax = 50;

N = 101; % number of samples for plotting



verbose = 0;
% --- CALCULATION and DISPLAY --- %



deg = 180/pi;
ar = ar/deg; % convert to radians for calculations
ao = ao/deg;
damin = damin/deg;
damax = damax/deg;
t = 1000*t; % convert to nm for calculations

C = cos(ar)/cos(ao);
Cs = sign(C);

disp(['wavelength: ' num2str(L) 'nm']);
disp(['thickness: ' num2str(t/1000) 'um']);
disp(['index of refraction: ' num2str(n) ', modulation: ' num2str(dn)]);
disp(['reference angle: ' num2str(ar*deg) ' degrees']);
disp(['object angle: ' num2str(ao*deg) ' degrees']);
disp(['angular range: ' ...
       num2str(damin*deg) ' to ' num2str(damax*deg) ' degrees']);
disp(['wavelength range: ' ...
       num2str(dLmin) ' to ' num2str(dLmax) ' nm']);
disp(['number of samples for plotting: ' num2str(N)]);
if Cs > 0
    disp(['transmission (C = ' num2str(C) ')']);
else
    disp(['reflection (C = ' num2str(C) ')']);
end

% q = 8*pi*t*sin(ar)^2/L;
% if q <= 10
%     disp(['The hologram is thin, Klein Parameter  =  ' q]);
% end


% check validity of theory
v =  pi*dn*t/(2*L*n*sqrt(Cs*cos(ar)*cos(ao)));
if Cs > 0
    valid = v*abs(tan(v));
else
    valid = v*abs(tanh(v));
end
disp(['validity: ' num2str(valid) ' should be << ' num2str(2*pi*t*n*cos(ao)/L)]);


% max diffraction efficiency
if Cs > 0
    effmax = sin(v)^2;
else
    effmax = tanh(v)^2;
end
disp(['max efficiency: ' num2str(effmax*100) '%']);


dL = linspace(dLmin, dLmax, N);
da = linspace(damin, damax, N);

xiL = -2*pi*n*t*dL*sin((ar-ao)/2)^2/(cos(ao)*L^2);
xia = pi*n*t*sin(ar-ao)*da/(L*cos(ao));


if Cs > 0 % (transmission)
    effL = sin(sqrt(xiL.^2 + v^2)).^2./(1 + xiL.^2/v^2);
    effa = sin(sqrt(xia.^2 + v^2)).^2./(1 + xia.^2/v^2);
else % (reflection)
    effL = 1./( 1 + (1 - xiL.^2/v^2)./sinh(sqrt(v^2 - xiL.^2)).^2 );
    effa = 1./( 1 + (1 - xia.^2/v^2)./sinh(sqrt(v^2 - xia.^2)).^2 );
end
etamaxL = max(effL);
etamaxa = max(effa);

if verbose
    disp('        wavelength dev (nm)    DE (%)      norm DE (%)');
    disp( cat(2,(1:N)', dL', 100*effL', 100*effL'/etamaxL) );
    disp('        angle dev (deg)    DE (%)      norm DE (%)');
    disp( cat(2,(1:N)', da'*deg, 100*effa', 100*effa'/etamaxa) );
end

figure(1);
plot(dL,effL*100,'.-');
title('Diffraction Efficiency vs \Delta\lambda');
xlabel('\Delta\lambda (nm)'); ylabel('diffraction efficiency (%)');

figure(2);
plot(dL,effL/etamaxL*100,'.-');
title('Normalized Diffraction Efficiency vs \Delta\lambda');
xlabel('\Delta\lambda (nm)'); ylabel('normalized diffraction efficiency (%)');

figure(3);
plot(da*deg,effa*100,'.-');
title('Diffraction Efficiency vs \Delta\theta');
xlabel('\Delta\theta (\circ)'); ylabel('diffraction efficiency (%)');

figure(4);
plot(da*deg,effa/etamaxa*100,'.-');
title('Normalized Diffraction Efficiency vs \Delta\theta');
xlabel('\Delta\theta (\circ)'); ylabel('normalized diffraction efficiency (%)');





