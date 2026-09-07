# Author Summary

Researchers use single-cell sequencing to build computer models that
predict whether a patient has a disease. These models often appear highly accurate.
We asked a simple question: are they detecting the disease itself, or are they
detecting how the samples happened to be collected and processed?

We studied two lupus datasets. In both, information about study
logistics alone, including when a sample was collected, which processing group it
belonged to and its sequencing quality, predicted disease without using any
biological measurement. When study design and disease are closely linked, an
analysis of that dataset alone cannot separate the two explanations.

We found that two standard remedies disagree sharply. Statistically subtracting
design effects removed far more apparent accuracy than comparing patients who were
processed under matched conditions. We verified that this disagreement remained
with a nonlinear adjustment method. A separate internal control showed why simply
detecting processing information in a model does not prove that the model used it
to predict disease.

We propose five checks that researchers can run before claiming that a
patient-level model has detected biology. Our code, data, software environments and
automated checks accompany the manuscript so others can apply the protocol.
