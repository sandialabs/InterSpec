# Tutorials

## An overview of using InterSpec for Analysis
See [brief_ana_overview_InterSpec_Oct2018.pdf](brief_analysis_intro/brief_ana_overview_InterSpec_Oct2018.pdf).
Data files:
* [Example1.pcf](brief_analysis_intro/spectra/Example1.pcf)
* [Example1_background.pcf](brief_analysis_intro/spectra/Example1_background.pcf)
* [Example2d_worked.n42](brief_analysis_intro/spectra/Example2d_worked.n42)
* [Example3_worked.n42](brief_analysis_intro/spectra/Example3_worked.n42)

## Making a Detector Response Function
Tutorial on making a detector response function in InterSpec: [detector_characterization_brief_20190619.pdf](make_drf/detector_characterization_brief_20190619.pdf). 
Further information (also available within InterSpec in its help system) is in [make_drf_help_20190619.pdf](make_drf/make_drf_help_20190619.pdf)
* Example data for HPGe (see [source_info.txt](make_drf/cal_data_HPGe/source_info.txt) for source and detector information):
  * [drf_cal_HPGe_Am241.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_Am241.pcf), [drf_cal_HPGe_Cd109.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_Cd109.pcf), [drf_cal_HPGe_Co60.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_Co60.pcf), [drf_cal_HPGe_Na22.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_Na22.pcf), [drf_cal_HPGe_Ba133.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_Ba133.pcf), [drf_cal_HPGe_Co57.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_Co57.pcf), [drf_cal_HPGe_Cs137.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_Cs137.pcf), [drf_cal_HPGe_U232.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_U232.pcf), [drf_cal_HPGe_background.pcf](make_drf/cal_data_HPGe/drf_cal_HPGe_background.pcf)
    * [source_info.txt](make_drf/cal_data_HPGe/source_info.txt)
* Example data for NaI:
  * [drf_cal_NaI_am241.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_am241.pcf), [drf_cal_NaI_cd109.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_cd109.pcf), [drf_cal_NaI_co60.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_co60.pcf), [drf_cal_NaI_na22.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_na22.pcf), [drf_cal_NaI_ba133.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_ba133.pcf), [drf_cal_NaI_co57.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_co57.pcf), [drf_cal_NaI_cs137.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_cs137.pcf), [drf_cal_NaI_u232.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_u232.pcf), [drf_cal_NaI_y88.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_y88.pcf), [drf_cal_NaI_22Na_18F.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_22Na_18F.pcf), [drf_cal_NaI_background.pcf](make_drf/cal_data_NaI_3x3/drf_cal_NaI_background.pcf)
    * [source_info.txt](make_drf/cal_data_NaI_3x3/source_info.txt)



## Example Problems
* Determine nuclide, and distance a source is inside a cargo container in a high scatter environment
  * Problem Setup: [cargo_container_question.pdf](example_problems/one_over_r2/problem_1/cargo_container/cargo_container_question.pdf)
  * Problem Solution: [cargo_container_solution.pdf](example_problems/one_over_r2/problem_1/cargo_container/cargo_container_solution.pdf)
  * Spectrum file with data: [cargo_container_3x3NaI.n42](example_problems/one_over_r2/problem_1/cargo_container/cargo_container_3x3NaI.n42)
* Determine nuclide, activity, and depth of a buried source from a single HPGe measurement
  * Problem Setup: [buried_source_question.pdf](example_problems/determine_activity_shielding/problem_1/buried_source_question.pdf)
  * Problem Solution: [buried_source_solution.pdf](example_problems/determine_activity_shielding/problem_1/buried_source_solution.pdf)
  * Spectrum file with data: [buried_source_40%_HPGe.n42](example_problems/determine_activity_shielding/problem_1/buried_source_40%_HPGe.n42)
* Determine Uranium enrichment and mass
  * Problem Setup: [uranium_problem_one_question.pdf](example_problems/uranium_enrichment/problem_1/uranium_problem_one_question.pdf)
  * Problem Solution: [uranium_problem_one_solution.pdf](example_problems/uranium_enrichment/problem_1/uranium_problem_one_solution.pdf)
  * Spectrum file with data: [uranium_40%_HPGe_15cm.n42](example_problems/uranium_enrichment/problem_1/uranium_40%_HPGe_15cm.n42)

