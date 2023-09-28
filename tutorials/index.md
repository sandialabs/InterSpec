# Tutorials

## Relative Efficiency Analysis 
Relative efficiency analyis allows you to determine enrichments and relative activities of nuclides in a spectrum, without knowing the detector response function, shielding, geometry, or distances.  This presentation provides an overview of one of InterSpecs relative efficiency tools.
* [20220922_InterSpec_RelEff_Peaks_SAND2022_15323TR.pdf](rel_eff_peaks/20220922_InterSpec_RelEff_Peaks_SAND2022_15323TR.pdf).



## InterSpec familiarization Apr 2022
* Part 1: Initial Familiarization (1 hour)
  * Presentation: [InterSpec_familiarization.pdf](familiarization_Apr2022/20220420_InterSpec_familiarization.pdf).
  * Example spectrum file: [example_problem_1.n42](familiarization_Apr2022/example_problem_1.n42).
  * DRFs of commonly encountered detectors (for InterSpec v1.0.9): [common_drfs.tsv](familiarization_Apr2022/common_drfs.tsv).
* Part 2: 3 More In-depth usage (3 hours)
  * Presentation: [InterSpec_more_in_depth_use.pdf](familiarization_Apr2022/20220427_InterSpec_more_in_depth_use.pdf).
  * Example spectra files:
    * Trinitite: Spectra provided by Dave Mercer of Los Alamos National Laboratory - they are sample "B" in the paper at https://arxiv.org/abs/2103.06240
      * [trinitite_sample_b.n42](familiarization_Apr2022/trinitite_sample_b.n42)
      * [trinitite_sample_b_background.n42](familiarization_Apr2022/trinitite_sample_b_background.n42)
      * My results fitting peaks: [Trinitite_Sample_B_peaks_fit.n42](familiarization_Apr2022/Trinitite_Sample_B_peaks_fit.n42)
    * Energy calibration practice: [th232_energy_cal_practice.n42](familiarization_Apr2022/th232_energy_cal_practice.n42)
    * Uranium enrichment [uranium_40%_HPGe_15cm_peaks_fit.n42](familiarization_Apr2022/uranium_40%_HPGe_15cm_peaks_fit.n42)
    * Surface contamination: [Eu152_surface_contamination_uDetective.n42](familiarization_Apr2022/Eu152_surface_contamination_uDetective.n42)


## Soil Contamination and Trace Sources
* [20211117_IAEA_HPGe_InterSpec_SAND2021-14557TR.pdf](contamination/20211117_IAEA_HPGe_InterSpec_SAND2021-14557TR.pdf).  
  The content in this PDF was presented as part of a IAEA webinar that can be watched [here](https://elearning.iaea.org/m2/course/view.php?id=1224) (requires IAEA NUCLEUS account).

## Introduction To Spectroscopy with InterSpec
* [20210315_spectroscopy_intro_InterSpec.pptx](spec_intro_March2021/20210315_spectroscopy_intro_InterSpec.pptx).
* [spectrum_files.zip](spec_intro_March2021/spectrum_files.zip).
* For a video of this presentation, see the [Nuclear Science and Security Consortium website](https://nssc.berkeley.edu/events/nssc-virtual-learning-series/) or on [YouTube](https://www.youtube.com/watch?v=xrwRYhVTC7Y).

## Overview Course Sep 2020
See:
* [20200915_InterSpec_intro_part1.pdf](intro_course_Sep2020/20200915_InterSpec_intro_part1.pdf).
* [20200915_InterSpec_intro_part2.pdf](intro_course_Sep2020/20200915_InterSpec_intro_part2.pdf).

Data files:
* All files [spectrum_files.zip](intro_course_Sep2020/spectrum_files.zip)


## An overview of using InterSpec for Analysis (Oct 2018)
See [brief_ana_overview_InterSpec_Oct2018.pdf](brief_analysis_intro/brief_ana_overview_InterSpec_Oct2018.pdf).

Data files:
* [Example1.pcf](brief_analysis_intro/spectra/Example1.pcf)
* [Example1_background.pcf](brief_analysis_intro/spectra/Example1_background.pcf)
* [Example2d_worked.n42](brief_analysis_intro/spectra/Example2d_worked.n42)
* [Example3_worked.n42](brief_analysis_intro/spectra/Example3_worked.n42)

## Making a Detector Response Function
Tutorial on making a detector response function in InterSpec: [detector_characterization_brief_20190619.pdf](make_drf/detector_characterization_brief_20190619.pdf). 

Further information (also available within InterSpec in its help system) is in [make_drf_help_20190619.pdf](make_drf/make_drf_help_20190619.pdf)
* Example data for HPGe: [drf_cal_data_HPGe.zip](make_drf/cal_data_HPGe/drf_cal_data_HPGe.zip).  See [source_info.txt](make_drf/cal_data_HPGe/source_info.txt) for source and detector information.

* Example data for NaI: [drf_cal_data_3x3_NaI.zip](make_drf/cal_data_NaI_3x3/drf_cal_data_3x3_NaI.zip). See [source_info.txt](make_drf/cal_data_NaI_3x3/source_info.txt) for source and detector info.
  
  
## Example Problems
* Determine nuclide, and distance a source is inside a cargo container in a high scatter environment
  * Problem Setup: [cargo_container_question.pdf](example_problems/one_over_r2/problem_1/cargo_container_question.pdf)
  * Problem Solution: [cargo_container_solution.pdf](example_problems/one_over_r2/problem_1/cargo_container_solution.pdf)
  * Spectrum file with data: [cargo_container_3x3NaI.n42](example_problems/one_over_r2/problem_1/cargo_container_3x3NaI.n42)
* Determine nuclide, activity, and depth of a buried source from a single HPGe measurement
  * Problem Setup: [buried_source_question.pdf](example_problems/determine_activity_shielding/problem_1/buried_source_question.pdf)
  * Problem Solution: [buried_source_solution.pdf](example_problems/determine_activity_shielding/problem_1/buried_source_solution.pdf)
  * Spectrum file with data: [buried_source_40%_HPGe.n42](example_problems/determine_activity_shielding/problem_1/buried_source_40%_HPGe.n42)
* Determine Uranium enrichment and mass
  * Problem Setup: [uranium_problem_one_question.pdf](example_problems/uranium_enrichment/problem_1/uranium_problem_one_question.pdf)
  * Problem Solution: [uranium_problem_one_solution.pdf](example_problems/uranium_enrichment/problem_1/uranium_problem_one_solution.pdf)
  * Spectrum file with data: [uranium_40%_HPGe_15cm.n42](example_problems/uranium_enrichment/problem_1/uranium_40%_HPGe_15cm.n42)


# Other useful resources
- **FRMAC Gamma Spectroscopist Knowledge Guide.**: an freely available, extensive, and excellent, guide to gamma spectroscopy; available [here (osti.gov)](https://www.osti.gov/biblio/1763003-frmac-gamma-spectroscopist-knowledge-guide-revision) and [here (nnss.gov)](https://www.nnss.gov/docs/docs_FRMAC/_FRMAC_GammaSpec_KnowledgeGuide_2019-08_UUR.pdf)
- **High Resolution Gamma-Ray Spectrometry Analyses For Normal Operations and Radiological Incident Response**: available [here (epa.gov)](https://www.epa.gov/sites/default/files/2020-07/documents/guide_for_high_resolution_gamma_spectrometry_analyses_camera_ready.pdf)
- The **Passive Nondestructive Assay of Nuclear Materials**: freely available, extensive, and deep gamma spectroscopy and theory guide, available [here (lanl.gov)](https://www.lanl.gov/orgs/n/n1/FMTTD/neut_mc/pdfs/LA_UR_90_0732.pdf).  Its 2007 addendum, available [here (lanl.gov)](https://www.lanl.gov/org/ddste/aldgs/sst-training/_assets/docs/PANDA%202007%20Addendum/PANDA%202007%20Addendum.pdf), contains a plethora of additional information and analysis techniques.
- [IAEA Nuclear Security Detection Sciences and Technology Webinars](https://elearning.iaea.org/m2/course/index.php?categoryid=248), which include a couple presentation on using InterSpec, with some further videos [here](https://nucleus.iaea.org/sites/nuclear-instrumentation/Pages/Portable-detectors_videos.aspx) (requires IAEA NUCLEUS account to access).
- Draft specification for representing spectra in a QR code or URI: [spectrum_in_a_qr_code_uur_latest.pdf](references/spectrum_in_a_qr_code_uur_latest.pdf) (implemented in InterSpec as of v1.0.11), and a quick overview: [20230829_spectra_in_a_QR-code_SAND2023-08778O.pdf](references/20230829_spectra_in_a_QR-code_SAND2023-08778O.pdf).