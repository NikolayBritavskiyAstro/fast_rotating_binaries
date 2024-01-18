rule unzip_MESA_output:
     input:
        "src/data/data.zip"
     script:
         "src/scripts/unzip_MESA_output.py"

rule plot_save_mesa_comparison_panels_fin:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/age_vsurf_comparison_resolution_test.pdf"
        "src/tex/figures/HR_primary.pdf"
        "src/tex/figures/LOGS_lg_tsync_HD191495.pdf"
        "src/tex/figures/age_vsurf_comparison_panel.pdf"
    script:
        "src/scripts/plot_save_mesa_comparison_panels_fin.py"

rule plot_save_mesa_comparison_panels_acc_fin:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/age_vsurf_acceleration_panel.pdf"
        "src/tex/figures/age_omega_comparison_panel.pdf"
    script:
        "src/scripts/plot_save_mesa_comparison_panels_acc_fin.py"

rule plot_e_fin:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/age_e_comparison_panel.pdf"
    script:
        "src/scripts/plot_e_fin.py"

rule hr_plot_paper2023_fin:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/HR_REAL_ROT_Paper2023.pdf"
    script:
        "src/scripts/hr_plot_paper2023_fin.py"

rule plot_save_mesa_comparison_wind_fin:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/age_vsurf_comparison_onestar.pdf"
        "src/tex/figures/age_vsurf_comparison_onestar_resolution_vrot.pdf"
        "src/tex/figures/age_vsurf_comparison_onestar_resolution.pdf"
    script:
        "src/scripts/plot_save_mesa_comparison_wind_fin.py"

rule all_hr_fin:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/age_vsurf_hr_all.pdf"
    script:
        "src/scripts/all_hr_fin.py"

rule plot_save_mesa_ce:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/r_alpha_HD25631.pdf"
        "src/tex/figures/r_alpha_HD191495.pdf"
        "src/tex/figures/r_alpha_HD46485.pdf"
    script:
        "src/scripts/plot_save_mesa_ce.py"

rule plot_save_mesa_individual_fin1p25:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/p_q_1p25days.pdf"
        "src/tex/figures/p_q_3days.pdf"
        "src/tex/figures/p_q_5days.pdf"
        "src/tex/figures/p_q_10days.pdf"
        "src/tex/figures/p_q_50days.pdf"
        "src/tex/figures/p_q_70days.pdf"
    script:
        "src/scripts/plot_save_mesa_individual_fin1p25.py"

rule plot_save_mesa_eta:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/eta_omega.pdf"
    script:
        "src/scripts/plot_save_mesa_eta.py"

rule plot_min_a_in:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/a_min_RL_24_1_7.pdf"
        "src/tex/figures/a_min_RL_15_1p5_3p5.pdf"
        "src/tex/figures/a_min_RL_7_1_5.pdf"
    script:
        "src/scripts/plot_min_a_in.py"

rule plot_P_unstable:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/P_unstable.pdf"
    script:
        "src/scripts/plot_P_unstable.py"

rule plot_save_mesa_individual_pos_res:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/hr_panel_res.pdf"
    script:
        "src/scripts/plot_save_mesa_individual_pos_res.py"

rule plot_save_mesa_pos_eta_res:
    input:
        "src/data/data.zip"
    output:
        "src/tex/figures/eta_omega_res.pdf"
    script:
        "src/scripts/plot_save_mesa_pos_eta_res.py"

