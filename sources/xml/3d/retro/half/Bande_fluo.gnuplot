set term epslatex standalone color font 8

set xlabel "incidence elevation (in degrees)"
set ylabel "BRDF"

# output ABC fits
set output "gray_retro_abc.tex"
plot "../papers/retro/mesures/original/Bande_fluo_grise/3d/633nm/Fichiers\ definitifs/densityHelmholtz/Bande_grise_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Gray cloth data", "./results/3d/retro/half/Bande_fluo_abc_retro.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro ABC fit", "./results/3d/retro/half/Bande_fluo_abc_back.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back ABC fit"

# output Beckmann fits
set output "gray_retro_beck.tex"
plot "../papers/retro/mesures/original/Bande_fluo_grise/3d/633nm/Fichiers\ definitifs/densityHelmholtz/Bande_grise_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Gray cloth data", "./results/3d/retro/half/Bande_fluo_beck_retro.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro Beckmann fit", "./results/3d/retro/half/Bande_fluo_beck_back.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back Beckmann fit"

# output Blinn fits
set output "gray_retro_blinn.tex"
plot "../papers/retro/mesures/original/Bande_fluo_grise/3d/633nm/Fichiers\ definitifs/densityHelmholtz/Bande_grise_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Gray cloth data", "./results/3d/retro/half/Bande_fluo_blinn_retro.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro Blinn fit", "./results/3d/retro/half/Bande_fluo_blinn_back.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back Blinn fit"

# output Lafotune fit
set output "gray_retro_laf.tex"
plot "../papers/retro/mesures/original/Bande_fluo_grise/3d/633nm/Fichiers\ definitifs/densityHelmholtz/Bande_grise_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Gray cloth data", "./results/3d/retro/half/Bande_fluo_laf.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Lafortune fit"

set ylabel "BRDF x cosine"
set size square

# output Lafotune fit
set output "gray_retro_rat.tex"
plot "../papers/retro/mesures/original/Bande_fluo_grise/3d/633nm/Fichiers\ definitifs/densityHelmholtz/Bande_grise_3D_dense__nbsgrid_162.alta" using (180/pi*$2):($3 > -0.01 && $3 < 0.05 ? $4 : 1/0) title "Gray cloth data", "./results/3d/retro/half/Bande_fluo_raf.dat" using (180/pi*$2):($3 > -0.01 && $3 < 0.05 ? $4 : 1/0) title "Rational interpolation"
