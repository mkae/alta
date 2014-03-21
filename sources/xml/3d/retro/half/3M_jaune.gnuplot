set term epslatex standalone color font 8

set xlabel "incidence elevation (in degrees)"
set ylabel "BRDF"
set key on inside left

# output ABC fits
set output "yellow_retro_abc.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_abc_retro.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro ABC fit", "./results/3d/retro/half/3M_jaune_abc_back.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back ABC fit"

# output Beckmann fits
set output "yellow_retro_beck.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_beck_retro.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro Beckmann fit", "./results/3d/retro/half/3M_jaune_beck_back.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back Beckmann fit"

# output Blinn fits
set output "yellow_retro_blinn.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_blinn_retro.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro Blinn fit", "./results/3d/retro/half/3M_jaune_blinn_back.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back Blinn fit"

# output Lafotune fit
set output "yellow_retro_laf.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_laf.dat" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Lafortune fit"

set ylabel "BRDF x cosine"

# output rational fits
set output "yellow_retro_rat.tex"
plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR__BRDF_min_retro_lobe_dense.alta" using (180/pi*$2):($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Yellow cloth data", "./results/3d/retro/half/3M_jaune_rat.dat" using (180/pi*$2):($3 > -0.01 && $3 < 0.05 ? $3 : 1/0) title "rational interpolation"
