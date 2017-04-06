# holography
Collection of files and scripts related to designing holograms and optical setups.

    beam_ratios.m:

Object and reference beams need to have a certain power ratio in the polymer. This script takes into account Fresnel reflections, beam density changes in different media, detector size, and polymer specs to give the power that the two beams should have in the air, as well as the total recording time.


    Kogelnik_Analysis.m:

Performs the analysis done by Kogelnik (1969 coupled-wave holography paper). Plots wavelength and angular selectivity for a specific hologram. Includes option to use take hologram specs from "Recording_Angles.m", but can be inputted manually also. Also includes option to use some linear approximates that Kogelnik makes which are actually not necessary; the plots can be made with the more exact forms.


    Recording_Angles.m:

Finds the optical geometry needed to record a hologram with one wavelength, when it is intended to be used at another wavelength, according to Bragg's Law. (You may not have the laser with the desired wavelength to record the hologram.) Designed to give angles in air, taking into account refractions into/out of the polymer and any index-matched substrate/glass/material around the polymer. Takes into account index dispersion, if the material properties are known.


    indexof.m:

Gives the index of refraction of a particular material at a given (or vector of) wavelenghts. Uses data from "www.refractiveindex.info" for most materials. Small corrections applied to the dispersion series for some materials in order to better fit data.


    thick_holo_huang.m:

Follows analysis done by Gibert and Huang in 1996 paper. Similar analysis to Kogelnik specifically for substrate guided-wave holograms. Practically, analysis is identical to Kogelnik, maybe even slightly worse. Probably best to use Kogelnik directly.
