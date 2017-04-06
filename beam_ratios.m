function beam_ratios(a1g,a1a,f1,a2g,a2a,f2,E,dA)
% recording beam ratio
% assuming we want power densities in glass to be equal

% a1g = 12.6576; % beam 1 angle in glass
% a1a = 19.3833; % beam 1 angle in air
% f1 = .10919; % fraction of power density lost to fresnel reflection
% 
% a2g = 29.6400; % beam 2 angle in glass
% a2a = 48.5076; % beam 2 angle in air
% f2 = .048735; % fraction of power density lost to fresnel reflection
% 
% E = 15; % energy density suggested by Covestro (mJ/cm2)
% dA = .708822; % detector area (cm2)


r1 = cos(a1g)/cos(a1a)/(1-f1);
r2 = cos(a2g)/cos(a2a)/(1-f2);
r = r1/r2;
E1 = E/2*r1*dA;
E2 = E/2*r2*dA;

p1a = [.1 .25 .5 .75 1 1.25 1.5]; % mW/detector area
p2a = p1a/r;
T1 = E1./p1a;
% T2 = E2./p2a; % should be same as T1



deg = 180/pi;

disp(['energy density suggested for polymer: ' num2str(E) ' mJ/cm2']);
disp(['detector area: ' num2str(dA) ' cm2']);
disp(' ');

disp('angles in glass (wrt normal):');
disp([a1g a2g]*deg);

disp('angles in air (wrt normal):');
disp([a1a a2a]*deg);

disp('fraction of power lost due to fresnel reflection (%):');
disp([f1 f2]*100);

disp('power density in air / in glass (beam 1, beam 2):');
disp([r1 r2]);

disp('for power densities in glass to be equal,');
disp('power density of beam 1 / power density of beam 2:');
disp(r);

disp('total energy per detector area (mJ) (beam 1, beam 2):');
disp([E1 E2]);

disp('example recording conditions:');
disp('beam 1 power per detector area (mW):');
disp(p1a);
disp('beam 2 power per detector area (mW):');
disp(p2a);
disp('associated recording time (seconds):');
disp(T1); % should be same at T2, no need to print both


end




