# holography
Collection of files and scripts related to designing holograms and optical setups.

-- beam_ratios.m --
Object and reference beams need to have a certain power ratio in the polymer. This script takes into account Fresnel reflections, beam density changes in different media, detector size, and polymer specs to give the power that the two beams should have in the air, as well as the total recording time.

-- Kogelnik_Analysis.m --
Performs the analysis done by Kogelnik (1969 coupled-wave holography paper). Plots wavelength and angular selectivity for a specific hologram. Includes option to use take hologram specs from "Recording_Angles.m", but can be inputted manually also. Also includes option to use some linear approximates that Kogelnik makes which are actually not necessary; the plots can be made with the more exact forms.
