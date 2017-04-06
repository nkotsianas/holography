% "Recording_Angles.m"
% Nick Kotsianas
% Version: 2016-05-31
%
% Dependencies:
%   (1) "indexof.m"
%
%
% This script finds the recording angles needed to record a hologram
%   at wavelength Lr when it is intended to be played back at Lp.
% Also optionally finds recording parameters (beam energy, recording time,
%   etc...) (in BETA!!!), and does Kogelnik analysis.


% --- USER INPUTS: --- %

Lp = 532; % playback wavelength (in nm)
    medlp = 'AIR'; % medium left of polymer (-90 to +90) at playback
    medrp = 'FUSEDSILICA'; % medium right of polymer (+90 to +180, -90 to -180)
Lr = 325; % recording wavelength (in nm)
    medlr = 'AIR'; % medium left of polymer (-90 to +90) at recording
    medrr = 'FUSEDSILICA'; % medium right of polymer (+90 to +180, -90 to -180)

recmedp = 'PTR-BK7'; % recording medium (ie holographic polymer) at playback
recmedr = recmedp; % recording medium (ie holographic polymer) at recording

% beam angles (in medium) (wrt normal) (for Lp nm beam):
Ramp = 0;
Samp = 43;



verbose = 1; % print out for angles in medium
kogelnik = 0; % do kogelnik plots and analysis
verboseair = 0; % print out for angles in air, beam ratio, fresnel
                % reflections, recording times, etc... (in BETA!)
% --------------------------------------------------- %
% --- everything below is calculation and display --- %
% --------------------------------------------------- %

% --- CALCULATION: --- %

% Use the function "indexof.m" to find all the refractive indices:
npp = indexof(recmedp,Lp); % index of polymer at playback
npr = indexof(recmedr,Lr); % index of polymer at recording
nmlp = indexof(medlp,Lp); % index of left medium at playback
nmlr = indexof(medlr,Lr); % index of left medium at recording
nmrp = indexof(medrp,Lp); % index of right medium at playback
nmrr = indexof(medrr,Lr); % index of right medium at recording

% - Define a few functions to make life easier - %
%   (and minimize external dependencies):

unit = @(V) V/norm(V); % Normalize a vector to its unit vector.

% Return the direction of a 2-element vector
ang = @(V) atan2(V(2),V(1));

% Snell's Law for vectors.
%   Requires V, N to be normalized and N points towards the incoming light.
%   Vector that is returned is also normalized.
snell = @(n1,n2,V,N) n1/n2*V + ...
    (-dot(V,N)*n1/n2 - sqrt(1 - (n1/n2)^2*(1 - (-dot(V,N))^2)))*N;


deg = 180/pi; % "A*deg" converts A from rad to deg

% -                                            - %

Ramp = Ramp/deg;
Samp = Samp/deg;



% Define the playback beam vectors in medium:
%   They can be unit length since it's not important to know their
%   magnitudes for any calculations.
Rmp = [cos(Ramp) sin(Ramp)];
Smp = [cos(Samp) sin(Samp)];


% Find the playback beam vectors in polymer:
%   Snell the reference beam into the polymer
%   (we assume the reference beam always comes from the left)
Rpp = snell(nmlp,npp,Rmp,[-1 0]);

%   Snell the signal (object) beam into the polymer:
if Smp(1) > 0 % in playback, the beam goes to the right (transmission)
    Spp = snell(nmrp,npp,Smp,[-1 0]);
else % in playback, the beam goes to the left (reflection)
    Spp = snell(nmlp,npp,Smp,[1 0]);
end

%   Once we have the directions, give them the appropriate magnitudes.
Rpp = (2*pi*npp/Lp) * Rpp;
Spp = (2*pi*npp/Lp) * Spp;


% Find the fringe vector:
%   Points along the direction perpendicular to the fringes.
%   Has magnitude 2*pi/d, where d is fringe spacing.
K = Rpp - Spp;

% Calculate fringe spacing (in units of Lp)
d = 2*pi/norm(K);

% Find the *parallel* fringe vector (points parallel to fringes):
%   Has magnitude sqrt(4*(2*pi*npp/Lp)^2 - (2*pi/d)^2).
K2p = Rpp + Spp;



% Calculate the *parallel* fringe vector for the recording beams:
%   Direction is same as K2p, but magnitude is calculated for recording.
K2r = sqrt(4*(2*pi*npr/Lr)^2 - (2*pi/d)^2) * unit(K2p);

% Calculate the recording beam vectors in polymer from K2r and K:
Rpr = (K2r + K)/2;
Spr = (K2r - K)/2;

% Snell the recording beams out of the polymer, into the surrounding medium
if Rpr(1) > 0 % Coming from the left
    Rmr = snell(npr,nmlr,unit(Rpr),[-1 0]);
else % Coming from the right
    Rmr = snell(npr,nmrr,unit(Rpr),[1 0]);
end

if Spr(1) > 0 % Coming from the left
    Smr = snell(npr,nmlr,unit(Spr),[-1 0]);
else % Coming from the right
    Smr = snell(npr,nmrr,unit(Spr),[1 0]);
end

% Once we have the directions, give them the appropriate magnitudes.
Rmr = (2*pi*nmlr/Lr) * Rmr;
Smr = (2*pi*nmlr/Lr) * Smr;



% % K-vector closure plots (well, attempts):
% figure(1);
% % plot([0 Rpp(1); 0 Spp(1); Spp(1) Spp(1)+K(1); 0 K2p(1)/2; Rpp(1)-Rpr(1) Rpp(1); Spp(1)-Spr(1) Spp(1)]', ...
% %      [0 Rpp(2); 0 Spp(2); Spp(2) Spp(2)+K(2); 0 K2p(2)/2; Rpp(2)-Rpr(2) Rpp(2); Spp(2)-Spr(2) Spp(2)]');
% % legend('Rpp', 'Spp', 'K', 'K2p/2', 'Rpr', 'Spr');
% plot([0 Rpp(1)],[0 Rpp(2)]);
%     text(Rpp(1)/2, Rpp(2)/2, 'Rpp', 'vert','t', 'horiz','l');
% hold on;
% plot([0 Spp(1)],[0 Spp(2)]);
%     text(Spp(1)/2, Spp(2)/2, 'Spp', 'vert','t', 'horiz','l');
% plot([Spp(1) Rpp(1)], [Spp(2) Rpp(2)],'r');
%     text(Spp(1)/2, Spp(2)/2, 'K', 'vert','t', 'horiz','l');
% rectangle('Position',[0-norm(Rpp),0-norm(Rpp),norm(Rpp)*2,norm(Rpp)*2], ...
%     'Curvature', [1,1]);
% rectangle('Position',[Rpp(1)-Rpr(1)-norm(Rpr),Rpp(2)-Rpr(2)-norm(Rpr),norm(Rpr)*2,norm(Rpr)*2], ...
%     'Curvature', [1,1]);
% hold off;
% axis equal;
% 
% figure(2);
% plotv([ Rpp' Spp' K' K2p' Rpr' Spr' ]);
% legend('Rpp', 'Spp', 'K', 'K2p', 'Rpr', 'Spr');
% axis equal;



% --- DISPLAY TO USER: --- %

if verbose
    disp(['Playback wavelength: ',num2str(Lp),'nm']);
    disp(['    Polymer index at playback: ',num2str(npp)]);
    disp(['    Left medium index at playback: ',num2str(nmlp),' (',medlp,')']);
    disp(['    Right medium index at playback: ',num2str(nmrp),' (',medrp,')']);
    disp(['Recording wavelength: ',num2str(Lr),'nm']);
    disp(['    Polymer index at recording: ',num2str(npr)]);
    disp(['    Left medium index at recording: ',num2str(nmlr),' (',medlr,')']);
    disp(['    Right medium index at recording: ',num2str(nmrr),' (',medrr,')']);
    disp(' ');

    disp([' ======== for playback (',num2str(Lp),'nm) beams ========']);
    
    disp(['beam angles (in medium) (wrt normal) (for ',num2str(Lp),'nm beam):']);
    disp([ Ramp Samp ]*deg);
    
    disp(['beam angles (in polymer) (wrt normal) (for ',num2str(Lp),'nm beam):']);
    disp([ ang(Rpp) ang(Spp) ]*deg);
    
    disp(['relative fringe angle (wrt ',num2str(Lp),'nm beams in polymer):']);
    disp(acosd(dot(unit(Rpp), unit(K2p))));
    
    disp('absolute fringe angle (wrt normal of polymer):');
    disp(ang(K2p)*deg);

    disp('fringe spacing (microns):');
    disp(d/1000);
    
    
    disp([' ======== for recording (',num2str(Lr),'nm) beams ========']);

    disp(['relative fringe angle (wrt ',num2str(Lr),'nm beams in polymer)']);
    disp(acosd(dot(unit(Rpr),unit(K2r))));

    disp(['beam angles (in polymer) (wrt normal) (for ',num2str(Lr),'nm beam):']);
    disp([ ang(Rpr) ang(Spr) ]*deg);

    disp(['beam angles (in medium) (wrt normal) (for ',num2str(Lr),'nm beam):']);
    disp([ ang(Rmr) ang(Smr) ]*deg);

    disp(['interbeam angle (in medium) (for ',num2str(Lr),'nm beams):']);
    disp(acosd(dot(unit(Rmr),unit(Smr))));
    
end






% --------------------------------------------------- %
% ---               K O G E L N I K               --- %
% --------------------------------------------------- %

if kogelnik
    disp(' -----  K O G E L N I K  ----- ');
    disp('Not yet!');
end























