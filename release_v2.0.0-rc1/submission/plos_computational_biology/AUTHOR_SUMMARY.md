# Author Summary

Single-cell RNA sequencing measures thousands of cells from each patient, but many
clinical studies ultimately predict one outcome per patient. A model can therefore
appear highly accurate because it captures disease biology, differences in how
samples were collected and processed, or both. We reconstructed sample-design
information for two public lupus cohorts and tested three common patient
representations: frozen Geneformer embeddings and two expression-based pseudobulk
summaries. Study-design variables alone predicted disease with high accuracy, and
all three molecular representations identified processing batches almost perfectly.
However, two standard attempts to diagnose this problem disagreed. Comparing
patients within overlapping batches largely preserved disease discrimination,
whereas mathematically removing design-associated components sometimes erased most
of it. The latter operation cannot distinguish technical bias from real biology when
study design and disease status are closely aligned. Independent transfer between
batches and cohorts exposed a more consistent loss of generalisation, and learned
pooling did not repair it. We propose four inexpensive checks that patient-level
single-cell classifiers should report: disease prediction from design alone, design
prediction from the representation, analysis within overlapping designs with
matched controls, and fully source-fitted external validation.

