function n = indexof(material, L)
% n = indexof(material, L)
% Please give L in nm.
% Available materials: 
%   AIR (sampled over 230-1690nm)
%   B270 (sampled over 435.83-656.27nm, rest assumed like BK7)
%   BAYFOL_CURED (sampled at 589.29nm, assumed like PMMA)
%   BAYFOL_UNREC (sampled at 589.29nm, assumed like PMMA)
%   BK7 (sampled over 300-2500nm)
%   FUSEDSILICA (sampled over 210-3710nm)
%   PMMA (sampled over 420-1620nm)
%   PTR-BK7 (sampled at 587.5 and 1085nm, rest assumed like BK7)
%   PTR-FS (sampled at 587.5 and 1085nm, rest assumed like FS)
% 
% Most data taken from "refractiveindex.info".
% Written by Nick Kotsianas
% Version 05-27-2016


L = L/1000; % convert L to mircon for calculation
material = upper(material);

switch material
    case 'PMMA'
        % sampled over 420-1620nm
        n = sqrt( 2.1778 + (6.1209e-3)*L.^2 - (1.5004e-3)*L.^4 + ...
            (2.3678e-2)*L.^-2 - (4.2137e-3)*L.^-4 + (7.3417e-4)*L.^-6 - ...
            (4.5042e-5)*L.^-8 ) + 1./(-950.26*L + 159.25);

    case 'BK7'
        % sampled over 300-2500nm
        n = sqrt( 1 + 1.03961212*L.^2./(L.^2 - 0.00600069867) + ...
            0.231792344*L.^2./(L.^2 - 0.0200179144) + ...
            1.01046945*L.^2./(L.^2 - 103.560653) );

    case 'BAYFOL_UNREC'
        % sampled at 589.29nm, assumed like PMMA
        n = indexof('PMMA',L*1000) - .005;

    case 'BAYFOL_CURED'
        % sampled at 589.29nm, assumed like PMMA
        n = indexof('PMMA',L*1000) + .002;
        
    case 'AIR'
        % sampled over 230-1690nm
        n = 1 + 0.05792105./(238.0185-L.^-2)+0.00167917./(57.362-L.^-2);
        
    case 'B270'
        % sampled over 435.83-656.27nm, rest assumed like BK7
        n = indexof('BK7',L*1000)*1.12 - ... 
            indexof('BK7',540.5)*(1.12-1) + .0065;
%         n = -.4991*L.^3 + .9908*L.^2 - .6923*L + 1.6889;
%         n = 10.001*L.^5 - 25.93*L.^4 + 26.189*L.^3 - 12.64*L.^2 + ...
%             2.7627*L + 1.3412;
    
    case 'FUSEDSILICA'
        % sampled over 210-3710nm
        n = sqrt( 1 + 0.6961663./(1 - 0.0684043^2*L.^-2) + ...
            0.4079426./(1 - 0.1162414^2*L.^-2)+...
            0.8974794./(1-9.896161^2*L.^-2) );
        
    case 'PTR-FS'
        % sampled at 587.5 and 1085nm, rest assumed like FS
        n = indexof('FUSEDSILICA',L*1000) + ...
            (1.4959 - indexof('FUSEDSILICA',587.5)) - ...
            (L*1000 - 587.5)/(1085-587.5)*1.14765e-04;
        
    case 'PTR-BK7'
        % sampled at 587.5 and 1085nm, rest assumed like BK7
        n = indexof('BK7',L*1000) + ...
            (1.4959 - indexof('BK7',587.5)) + ...
            (L*1000 - 587.5)/(1085-587.5)*1.24212661e-03;

end


end







