set term epslatex standalone color font 8

set xlabel "incidence elevation (in degrees)"
set ylabel "BRDF"

# output ABC fits
set output "orange_retro_abc.tex"
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_abc_retro.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "retro ABC fit" w points pointtype 12, "./results/3d/retro/half/Bande_orange_abc_back.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "back ABC fit"

# output Beckmann fits
set output "orange_retro_beck.tex"
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_beck_retro.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "retro Beckmann fit" w points pointtype 12, "./results/3d/retro/half/Bande_orange_beck_back.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "back Beckmann fit"

# output Blinn fits
set output "orange_retro_blinn.tex"
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_blinn_retro.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "retro Blinn fit" w points pointtype 12, "./results/3d/retro/half/Bande_orange_blinn_back.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "back Blinn fit"

# output Lafotune fit
set output "orange_retro_laf.tex"
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_laf.dat" using (180/pi*$2):($3 > -0.001 && $3 < 0.001 ? $4 : 1/0) title "Lafortune fit"

set ylabel "BRDF x cosine"

# output Rational fit
set term wxt
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D_dense__nbsgrid_162.alta" using (180/pi*$2):($2 > 0.0 && $3 > -0.01 && $3 < 0.05 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_rat_TK2TL.dat" using (180/pi*$2):($2 > 0.0 && $3 > -0.01 && $3 < 0.05 ? $4 : 1/0) title "Rational interpolation"

set term epslatex standalone color font 8
set output "orange_retro_rat.tex"
set yrange [0:GPVAL_Y_MAX]
replot
