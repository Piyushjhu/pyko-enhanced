# CHANGELOG pyKO example inputs

## 2026-07-20
* Fixed a unit-conversion bug in Test 18a/18b (`pos` was scaled by 1e4 instead of 1e6 m→μm) that distorted position axes in FSV plots.
* Added a `near_existing_fracture()` guard in `pyko.py` to stop runaway re-fracturing: once a zone spalls or an interface separates, neighboring zones within a radius are no longer re-fractured every time step.
* Added optional top-level `impact_velocity` YAML key that sets the mat1 flyer speed directly, instead of requiring `mat1.init.up0`; `material_database.py` now generates configs using this key.
* Test 18a: FSV plots are now aligned to shock-arrival time, with a zoomed view centered on arrival instead of a fixed time window.
* Test 18b: impact velocity is now read from `thickness_sweep.constant_velocity` (previously `mat1.init.up0`), and the script errors clearly if it's missing; removed unused minimum-`tstop` check.
* Test 18a/18b example configs (`test18.yml`, `test18b.yml`) switched target material from Brass/Aluminum to Copper and updated velocity/thickness sweep ranges.
* `test17-spall-interface` scripts now resolve paths relative to the script location instead of a hard-coded absolute path, so they run from any working directory.
* Added `pko-test17-xt-animation.py`, an x-t diagram animation script for Test 17.

## 2023-07-05
* Added Test8 notebook for a lab experiment on ice and two tabular water EOS that were not generated with ANEOS: the 5PHASE EOS and AQUA EOS. 
  See those eos directories for how to convert a table into a format compatible with pyKO.

## 2023-07-03
* Added Test13 notebooks to compare Tillotson, Mie Grueneisen, and tabular EOS. Examples for Al and Fe. Now including temperature calculations.
* Added Test14 notebook to compare Hosono, iSALE, and Asphaug implementations of Tillotson. Cross check of calculated sound speeds.

## 2023-06-28
* Added x-t diagram to test1b jupyter notebook

## v0.6.0 2023-06-27

Example input files are included with Jupyter Notebooks, for training and visualization.

* Test 1: Two Mie-Grueneisen plates planar impact, comparison to fortran KO
* Test 1b: Three MGR plates
* Test 2: Sod shock tube test for ideal gases, comparison to analytic solution and fortran KO
* Test 3: EOS table ideal gas Sod shock tube test, comparison to analytic solutoin
* Test 4: EOS table planar plate impact, comparison to Mie-Grueneisen fortran KO, includes movie visualization examples
* Test 6: MGR plates with gaps, example of closing gaps
* Test 6b: SES plates with gaps, example of closing gaps
* Test 9: Spall test with MGR plates
* Test 9b: Spall test with TIL and MGR plates
* Test 12: Gravity Al target and impactor
* Test 12b: Gravity Forsterite target and impactor
