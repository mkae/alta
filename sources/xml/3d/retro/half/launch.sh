#! /bin/sh

# Launch fits
./scripts/xml_cmd.py ./xml/3d/retro/half/3M_jaune_abc.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/3M_jaune_laf.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/3M_jaune_beck.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/3M_jaune_blinn.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/3M_jaune_rat.xml

./scripts/xml_cmd.py ./xml/3d/retro/half/Bande_fluo_abc.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/Bande_fluo_laf.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/Bande_fluo_beck.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/Bande_fluo_blinn.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/3M_jaune_blinn.xml

./scripts/xml_cmd.py ./xml/3d/retro/half/Bande_orange_abc.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/Bande_orange_laf.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/Bande_orange_beck.xml
./scripts/xml_cmd.py ./xml/3d/retro/half/Bande_orange_blinn.xml

# Generate latex files with figures
gnuplot ./xml/3d/retro/half/3M_jaune.gnuplot
gnuplot ./xml/3d/retro/half/Bande_fluo.gnuplot
gnuplot ./xml/3d/retro/half/Bande_orange.gnuplot

# Compile figures to pdf files
pdflatex yellow_retro_abc.tex
pdflatex yellow_retro_laf.tex
pdflatex yellow_retro_beck.tex
pdflatex yellow_retro_blinn.tex
pdflatex yellow_retro_rat.tex

pdflatex gray_retro_abc.tex
pdflatex gray_retro_laf.tex
pdflatex gray_retro_beck.tex
pdflatex gray_retro_blinn.tex

pdflatex orange_retro_abc.tex
pdflatex orange_retro_laf.tex
pdflatex orange_retro_beck.tex
pdflatex orange_retro_blinn.tex

rm *-inc.eps *.log *.tex *.aux *-to.pdf
