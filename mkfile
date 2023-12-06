OCTL=octave-cli
IMGSTATS=imgstats

IMGSTATSBIN=/home/me/go/bin/imgstats
VORFILTBIN=/home/me/go/bin/voronoi_filters

EXP_1_2_3_IMGS=`{ls ./images/exp1/*.png}
EXP11A_IMGS=`{ls ./images/exp11a/*.png}
EXP11B_IMGS=`{ls ./images/exp11b/*.png}
EXP5_IMGS=`{ls ./images/exp5/*.png}
RMC_RSD_IMGS=`{ls ./images/exp5_rmc_rsd_comp/*.png}

IMG_STATS=`{ls ./data/*_img_stats.csv}
OBS_STATS=`{ls ./data/*_obs_data.csv}

FIGS=./figures/refl_illum_coords.pdf ./figures/all_exps_mean_lab.pdf ./figures/white_wall_results_exp11a.pdf ./figures/white_wall_results_exp11b.pdf ./figures/cer_and_sd_results.pdf ./figures/rmc_rsd_comp_EXTREMMEEE.pdf ./figures/rmc_rsd_from_simulation.pdf ./figures/rmc_rsd_comp_EXTREMMEEE_dark_pxs_excl.pdf

all:V: figs

view:V:
	mupdf ./paper/paper_formatted.pdf

paper:V: paper_formatted.pdf $FIGS

paper_formatted.pdf:
	cd ./paper
	pdflatex paper_formatted.tex
	bibtex paper_formatted
	pdflatex paper_formatted.tex
	pdflatex paper_formatted.tex

figs:V: $FIGS

./figures/refl_illum_coords.pdf: ./analysis/plotIllumAndObjectColors.R
	cd ./analysis
	Rscript plotIllumAndObjectColors.R

./figures/all_exps_mean_lab.pdf ./figures/white_wall_results_exp11a.pdf ./figures/white_wall_results_exp11b.pdf: $OBS_STATS ./analysis/general_analysis_for_paper.R
	cd ./analysis
	Rscript general_analysis_for_paper.R

./figures/cer_and_sd_results.pdf: $OBS_STATS ./analysis/exp3.R
	cd ./analysis
	Rscript exp3.R

./figures/rmc_rsd_comp.pdf ./figures/rmc_rsd_from_simulation.pdf: $OBS_STATS $IMG_STATS ./analysis/exp5.R
	cd ./analysis
	Rscript exp5.R

./figures/rmc_rsd_comp_EXTREMMEEE_dark_pxs_excl.pdf: $OBS_STATS $IMG_STATS ./analysis/exp5_dark_pxs_excl.R
	cd ./analysis
	Rscript exp5_dark_pxs_excl.R

./data/%_obs_data.csv: ./data/%_img_stats.csv ./analysis/%.m
	cd ./analysis
	$OCTL $stem.m

./data/%_img_stats_dark_pxs_excl.csv: $IMGSTATSBIN $VORFILTBIN $EXP_1_2_3_IMGS $EXP11A_IMGS $EXP11B_IMGS $EXP5_IMGS $RMC_RSD_IMGS
	cd ./images/$stem/
	$IMGSTATS '-monName='$stem '-spectraFile='../../calibration/$stem.spectra.csv '-chromaFile='../../calibration/$stem.chroma.csv '-gammaFile='../../calibration/$stem.gamma.csv '-maskDir=../masks/' '-darkThresh=0.5' *.png
	mv $stem'_img_stats_dark_pxs_excl.csv' ../../data

./data/%_img_stats.csv: $IMGSTATSBIN $VORFILTBIN $EXP_1_2_3_IMGS $EXP11A_IMGS $EXP11B_IMGS $EXP5_IMGS $RMC_RSD_IMGS
	cd ./images/$stem/
	$IMGSTATS '-monName='$stem '-spectraFile='../../calibration/$stem.spectra.csv '-chromaFile='../../calibration/$stem.chroma.csv '-gammaFile='../../calibration/$stem.gamma.csv '-maskDir=../masks/' *.png
	mv $stem'_img_stats.csv' ../../data

$IMGSTATSBIN:
	cd /home/me/go/src/github.com/rennis250/imgstats/
	go install

$VORFILTBIN:
	cd /home/me/go/src/github.com/rennis250/voronoi_filters/
	go install