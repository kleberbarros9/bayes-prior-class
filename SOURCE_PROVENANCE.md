# Source provenance

Repository supplied by the authors: https://github.com/kleberbarros9/bayes-prior-class.git
Inspected initial commit: 7bee27ffb9579e5cc689d9b1c30d440b6436ded9 (README only).
This release consolidates the manuscript validation and final Roever scripts and the completed full lung and P-spline implementations supplied by the authors.
The final dataset, dose and transformations are fixed. No candidate search is executed.
The common survival Stan source removes the superseded comparator while retaining the same final-model likelihoods and common-target AFT priors. Consequently fresh random-number trajectories can differ from archived chains even with matching seeds.
Reference CSV files document previous numerical estimates; they are not substituted for newly computed outputs.
