rule shift_stack_uranus:
    input:
        "src/data/reduced/2019oct28/"
    output:
        "src/data/results/urh_Mab_2019-10-28.fits","src/data/results/urk_Mab_2019-10-28.fits"
    script:
        "src/scripts/shift_stack_uranus.py"

rule mab_path_on_detector:
    input:
        "src/data/reduced/2019oct28/"
    output:
        "src/tex/figures/motion_on_detector.png"
    script:
        "src/scripts/mab_path_on_detector.py"

rule detection_pretty_picture:
    input:
        "src/data/results/urh_Mab_2019-10-28.fits","src/data/results/urk_Mab_2019-10-28.fits"
    output:
        "src/tex/figures/detection_images_Mab.png"
    script:
        "src/scripts/detection_pretty_picture.py"

rule run_perturbation_experiment:
    input:
        "src/data/reduced/2019oct28/"
    output:
        "src/data/perturbation_experiment/urh_Mab_2019-10-28_99.fits","src/data/perturbation_experiment/urk_Mab_2019-10-28_99.fits"
    script:
        "src/scripts/run_perturbation_experiment.py"

rule plot_perturbation_experiment:
    input:
        "src/data/perturbation_experiment/urh_Mab_2019-10-28_99.fits","src/data/perturbation_experiment/urk_Mab_2019-10-28_99.fits"
    output:
	            "src/tex/figures/random_stack_experiment_Mab_H.png","src/tex/figures/random_stack_experiment_Mab_K.png","src/tex/output/perturbation_percent_higher_K.txt",
    script:
        "src/scripts/plot_perturbation_experiment.py"

rule psf_wing_correction:
    input:
        "src/data/reduced/2019oct28/"
    output:
        "src/data/tables/wing_correction_H_0.csv", "src/data/tables/wing_correction_H_1.csv", "src/data/tables/wing_correction_H_2.csv", "src/data/tables/wing_correction_K_0.csv", "src/data/tables/wing_correction_K_1.csv", "src/data/tables/wing_correction_K_2.csv"
    script:
        "src/scripts/psf_wing_correction.py"

rule photometry:
    input:
        "src/data/results/urh_Mab_2019-10-28.fits","src/data/results/urk_Mab_2019-10-28.fits", "src/data/tables/wing_correction_H_0.csv", "src/data/tables/wing_correction_H_1.csv", "src/data/tables/wing_correction_H_2.csv", "src/data/tables/wing_correction_K_0.csv", "src/data/tables/wing_correction_K_1.csv", "src/data/tables/wing_correction_K_2.csv"
    output:
        "src/tex/figures/photometry_vs_region_Mab_urh.png","src/tex/figures/photometry_vs_region_Mab_urk.png","src/tex/output/Mab_urh_flux.txt","src/tex/output/Mab_urk_fluxerr.txt"
    script:
        "src/scripts/photometry.py"

rule compute_if:
    input:
        "src/static/h.csv","src/static/kp.csv"
    output:
        "src/tex/output/color_table.tex","src/tex/figures/reflectance_spectrum.png","src/tex/output/Mab_urh_intif.txt","src/tex/output/Mab_urk_intif.txt","src/tex/output/Mab_urh_intiferr.txt","src/tex/output/Mab_urk_intiferr.txt"
    script:
        "src/scripts/compute_if.py"