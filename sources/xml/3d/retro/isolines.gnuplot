set linetype 1 linecolor rgb "red" linewidth 1 pointtype 1
set linetype 2 linecolor rgb "green" linewidth 1 pointtype 4

set term epslatex standalone color font 8 clip

#plot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR_dense__nbsgrid_162.alta" u 2:($1 < 0.4 && $3 > 0 && $3 < 0.02 && $2 > 0) ? $4 : 1/0

set dgrid3d 128,128,4
set view map
set size ratio -1
set contour
set cntrparam level incremental 2, 2, 8
set style data lines
set xrange [0:90]
set xlabel "$\\theta \\cos(\\Delta\\phi)$ (in degrees)"
set yrange [-15:15]
set ylabel "$\\theta \\sin(\\Delta\\phi)$ (in degrees)"
unset surface
unset clabel
#set output "yellow_isoline_classical.tex"
#splot "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR_dense__nbsgrid_162.alta" u (180/pi*$2):(180/pi*$3):($1 < 0.4 && $2 > 0) ? $4 : 1/0 notitle lt 1,  "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR_dense__nbsgrid_162.alta" u (180/pi*$2):(180/pi*$3):($1 < 0.9 && $1 > 0.4 && $2 > 0) ? $4 : 1/0 notitle lt 1, "../papers/retro/mesures/original/3M_jaune/3d/633nm/Fichiers definitifs/densify_helmholtz/3M_jaune_3D+3DS+3DR_dense__nbsgrid_162.alta" u (180/pi*$2):(180/pi*$3):($1 > 0.9 && $2 > 0) ? $4 : 1/0 notitle lt 1

#set output "yellow_isoline_back.tex"
#splot "/tmp/yellow_back.alta" u (180/pi*$2 + 15):(180/pi*$3):($1 < 0.4) ? $4 : 1/0 notitle lt 1,  "/tmp/yellow_back.alta" u (180/pi*$2 + 30):(180/pi*$3):($1 < 0.9 && $1 > 0.4) ? $4 : 1/0 notitle lt 1, "/tmp/yellow_back.alta" u (180/pi*$2 + 60):(180/pi*$3):($1 > 0.9) ? $4 : 1/0 notitle lt 1

#set output "yellow_isoline_retro.tex"
#splot "/tmp/yellow_retro.alta" u (180/pi*$2 + 15):(180/pi*$3):($1 < 0.4) ? $4 : 1/0 notitle lt 1,  "/tmp/yellow_retro.alta" u (180/pi*$2 + 30):(180/pi*$3):($1 < 0.9 && $1 > 0.4) ? $4 : 1/0 notitle lt 1, "/tmp/yellow_retro.alta" u (180/pi*$2 + 60):(180/pi*$3):($1 > 0.9) ? $4 : 1/0 notitle lt 1

set output "yellow_isoline_comparison.tex"
splot  "/tmp/yellow_back.alta" u (180/pi*$2 + 15):(180/pi*$3):($1 < 0.4) ? $4 : 1/0 notitle lt 1,  "/tmp/yellow_back.alta" u (180/pi*$2 + 30):(180/pi*$3):($1 < 0.9 && $1 > 0.4) ? $4 : 1/0 t "" lt 1, "/tmp/yellow_back.alta" u (180/pi*$2 + 60):(180/pi*$3):($1 > 0.9) ? $4 : 1/0 notitle lt 1, "/tmp/yellow_retro.alta" u (180/pi*$2 + 15):(180/pi*$3):($1 < 0.4) ? $4 : 1/0 notitle lt 2,  "/tmp/yellow_retro.alta" u (180/pi*$2 + 30):(180/pi*$3):($1 < 0.9 && $1 > 0.4) ? $4 : 1/0 notitle lt 2, "/tmp/yellow_retro.alta" u (180/pi*$2 + 60):(180/pi*$3):($1 > 0.9) ? $4 : 1/0 notitle lt 2

# Mac OS there is no need for waiting.
#pause -1
