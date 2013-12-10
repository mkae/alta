set term epslatex standalone color font 8

# output ABC fits
set output "orange_retro_abc.tex"
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D__BRDF_min_retro_lobe_dense.alta" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_abc_retro.dat" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro ABC fit", "./results/3d/retro/half/Bande_orange_abc_back.dat" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back ABC fit"

# output Beckman fits
set output "orange_retro_beck.tex"
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D__BRDF_min_retro_lobe_dense.alta" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_beck_retro.dat" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro Beckman fit", "./results/3d/retro/half/Bande_orange_beck_back.dat" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back Beckman fit"

# output Blinn fits
set output "orange_retro_blinn.tex"
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D__BRDF_min_retro_lobe_dense.alta" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_blinn_retro.dat" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "retro Blinn fit", "./results/3d/retro/half/Bande_orange_blinn_back.dat" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "back Blinn fit"

# output Lafotune fit
set output "orange_retro_laf.tex"
plot "../papers/retro/mesures/original/Bande_orange/3d/633nm/Fichiers_definitifs/densify_helmholtz/Bande_orange_3D__BRDF_min_retro_lobe_dense.alta" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Orange cloth data", "./results/3d/retro/half/Bande_orange_laf.dat" using 2:($3 > 0.0 && $3 < 0.005 ? $4 : 1/0) title "Lafortune fit"
