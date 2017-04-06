% "Kogelnik_Analysis.m"
% Nick Kotsianas
% Version: 2016-05-31
%
% Dependencies:
%   (1) "indexof.m"
%
%
% This script finds the angluar and wavelength selectivity of a
%   holographic grating using Kogelnik's coupled wave theory.
% This script is designed to take output from the "Recording_Angles.m"
%   file, but the data can be input manually also.


% Kogelnik makes some unnecessary linear approximations of a few quantities
%   (the dephasing factor, the two obliquity factors, and the two
%   parameters ni and xi).
% Put value 1 to do as Kogelnik does, or put 0 to use the full form of
%   those quantities.
doaskogelnikdoes = 1;


L = Lp; % incident wavelength (nm)
% d = 0.192975793480083; % fringe spacing (micron)

phi = ang(K); % angle of grating vector (*not* the fringes!)
a = ang(Rpp); % angle of incidence of reference beam in the *polymer*

% recmed = 'PTR-BK7'; % the recording medium (see "indexof.m" for materials)
dn = 0.001; % index variation amplitude
t = 1000; % polymer thickness (microns)

damin = -.05; % lower angular deviation (degrees)
damax = .05; % upper angular deviation (degrees)

dLmin = -.5; % lower wavelength deviation (nm)
dLmax = .5; % upper wavelength deviation (nm)

N = 1001; % samples for plotting (typically an odd number)





debug = 0; % makes plots of a few more parameters
% ----------------------------------------------------- %
% deg = 180/pi; % defined in "Recording_Angles.m"
% phi = phi/deg; % already in radians from "ang" in "Recording_Angles.m"
% a = a/deg; % already in radians from "ang" in "Recording_Angles.m"
damin = damin/deg;
damax = damax/deg;
t = t*1000;
% d = d*1000; % already in nanometers from "Recording_Angles.m"

da = linspace(damin, damax, N);
dL = linspace(dLmin, dLmax, N);


if doaskogelnikdoes
    
    n = indexof(recmed, L);

    deph = @(vda,vdL) vda*2*pi/d*sin(phi-a) - vdL*pi/(n*d^2); % dephasing parameter
    Cr = cos(a); % obliquity factor for reference beam
    Cs = cos(a) - L*cos(phi)/(n*d); % obliquity factor for signal beam
    ni = pi*dn*t/( L*sqrt(Cr*abs(Cs)) );
    xi = @(vda,vdL) t/2*deph(vda,vdL)/abs(Cs);

    if Cs > 0 % transmission
        eta = @(da,dL) sin(sqrt(ni^2 + xi(da,dL).^2)).^2 ./ ...
            (1 + (xi(da,dL)/ni).^2);
        etamax = sin(ni)^2;
    else % reflection
        eta = @(da,dL) 1./( 1 + (1 - (xi(da,dL)/ni).^2) .* ...
            csch(sqrt(ni^2 - xi(da,dL).^2)).^2 );
        etamax = tanh(ni)^2;
    end
    
    
    disp(' -----  K O G E L N I K  ----- ');
    disp(['Using Kogelnik''s linear approximations of the' ...
          ' dephasing factor, Cr, Cs, ni, and xi.']);
    disp(' ');

    disp('Wavelength (nm)');
    disp(L);
    disp(['Index of refraction at ' num2str(L) 'nm (' recmed ')']);
    disp(n);
    disp('Index modulation');
    disp(dn);
    disp('Grating spacing (microns)');
    disp(d/1000);
    disp('Grating vector angle (degrees)');
    disp(phi*deg);
    disp('Angle of incidence of reference beam in the polymer (degrees)');
    disp(a*deg);
    disp('Polymer thickness (micron)');
    disp(t/1000);
    
    disp('Bragg line (degrees/10nm)');
    disp(10/(2*d*n*sin(phi-a))*deg);
    
    disp('Max efficiency (%)');
    disp(etamax*100);
    
    
    
    figure(1);
    plot(da*deg, 100*eta(da,0));
    title('D.E. vs \Delta\theta');
    xlabel('\Delta\theta (\circ)'); ylabel('\eta (%)');
    ylim([0 100]);
    figure(2);
    plot(dL, 100*eta(0,dL));
    title('D.E. vs \Delta\lambda');
    xlabel('\Delta\lambda (nm)'); ylabel('\eta (%)');
    ylim([0 100]);
    
    disp(['ni = ' num2str(ni)]);
    disp(['Cr = ' num2str(Cr)]);
    disp(['Cs = ' num2str(Cs)]);
    
    % mostly for debugging:
    if debug
        figure(3);
        plot(da*deg, deph(da,0));
        title('dephasing factor vs \Delta\theta');
        xlabel('\Delta\theta (\circ)'); ylabel('\psi (nm^{-1})');
        figure(4);
        plot(dL, deph(0,dL));
        title('dephasing factor vs \Delta\lambda');
        xlabel('\Delta\lambda (nm)'); ylabel('\psi (nm^{-1})');

        figure(5);
        plot(da*deg, xi(da,0));
        title('\xi vs \Delta\theta');
        xlabel('\Delta\theta (\circ)');
        figure(6);
        plot(dL, xi(0,dL));
        title('\xi vs \Delta\lambda');
        xlabel('\Delta\lambda (nm)');
    end
    
    
    
    
else
    
    n = @(vL) indexof(recmed, vL);

    deph = @(vda,vdL) 2*pi/d*cos(a+vda-phi) - pi*(L+vdL)./(n(L+vdL)*d^2); % dephasing parameter
    Cr = @(vda) cos(a+vda); % obliquity factor for reference beam
    Cs = @(vda,vdL) cos(a+vda) - (L+vdL)*cos(phi)./(n(L+vdL)*d); % obliquity factor for signal beam
    ni = @(vda,vdL) pi*dn*t./( (L+vdL).*sqrt(Cr(vda).*abs(Cs(vda,vdL))) );
    xi = @(vda,vdL) t/2*deph(vda,vdL)./abs(Cs(vda,vdL));

    if Cs(0,0) > 0 % transmission
        eta = @(da,dL) sin(sqrt(ni(da,dL).^2 + xi(da,dL).^2)).^2 ./ ...
            (1 + (xi(da,dL)./ni(da,dL)).^2);
        etamax = sin(ni(0,0))^2;
    else % reflection
        eta = @(da,dL) 1./( 1 + (1 - (xi(da,dL)./ni(da,dL)).^2) .* ...
            csch(sqrt(ni(da,dL).^2 - xi(da,dL).^2)).^2 );
        etamax = tanh(ni(0,0))^2;
    end
    
    
    disp(' -----  K O G E L N I K  ----- ');
    disp(['Using full form of the dephasing factor, Cr, Cs, ni, and xi' ...
          ' (not Kogelnik''s linear approximations).']);
    disp(' ');

    disp('Wavelength (nm)');
    disp(L);
    disp(['Index of refraction at ' num2str(L) 'nm (' recmed ')']);
    disp(n(L));
    disp('Index modulation');
    disp(dn);
    disp('Grating spacing (microns)');
    disp(d/1000);
    disp('Grating vector angle (degrees)');
    disp(phi*deg);
    disp('Angle of incidence of reference beam in the polymer (degrees)');
    disp(a*deg);
    disp('Polymer thickness (micron)');
    disp(t/1000);
    
    disp('Bragg line (degrees/10nm)');
    disp(10/(2*d*n(L)*sin(phi-a))*deg);
    
    disp('Max efficiency (%)');
    disp(etamax*100);
    
    
    
    figure(1);
    plot(da*deg, 100*eta(da,0));
    title('D.E. vs \Delta\theta');
    xlabel('\Delta\theta (\circ)'); ylabel('\eta (%)');
    ylim([0 100]);
    figure(2);
    plot(dL, 100*eta(0,dL));
    title('D.E. vs \Delta\lambda');
    xlabel('\Delta\lambda (nm)'); ylabel('\eta (%)');
    ylim([0 100]);
    
    
    % mostly for debugging:
    if debug
        figure(3);
        plot(da*deg, deph(da,0));
        title('dephasing factor vs \Delta\theta');
        xlabel('\Delta\theta (\circ)'); ylabel('\psi (nm^{-1})');
        figure(4);
        plot(dL, deph(0,dL));
        title('dephasing factor vs \Delta\lambda');
        xlabel('\Delta\lambda (nm)'); ylabel('\psi (nm^{-1})');

        figure(5);
        plot(da*deg, ni(da,0), da*deg, xi(da,0));
        title('\nu, \xi vs \Delta\theta');
        xlabel('\Delta\theta (\circ)'); legend('\nu_\theta', '\xi_\theta');
        figure(6);
        plot(dL, ni(0,dL), dL, xi(0,dL));
        title('\nu, \xi vs \Delta\lambda');
        xlabel('\Delta\lambda (nm)'); legend('\nu_\lambda', '\xi_\lambda');

        figure(7);
        plot(da*deg, Cr(da), da*deg, Cs(da,0));
        title('Cr, Cs vs \Delta\theta');
        xlabel('\Delta\theta (\circ)'); legend('Cr_\theta', 'Cs_\theta');
        figure(8);
        plot(dL, ones(size(dL))*Cr(0), dL, Cs(0,dL));
        title('Cr, Cs vs \Delta\lambda');
        xlabel('\Delta\lambda (nm)'); legend('Cr_\lambda', 'Cs_\lambda');
    end

end


% [DA DL] = meshgrid(da,dL);
% figure();
% surf(DA*deg,DL, eta(DA,DL));


