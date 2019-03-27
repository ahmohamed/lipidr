# Internal script to update the dataset, in case of future changes.
# Should be called once before each package release.
library(lipidr)
datadir = system.file("extdata", package="lipidr")
filelist = list.files(datadir, "data.csv", full.names = TRUE)
d = read_skyline(filelist)
clinical_file = system.file("extdata", "clin.csv", package="lipidr")
d = add_sample_annotation(d, clinical_file)
d_summarized = summarize_transitions(d, method = "average")
data_normalized = normalize_pqn(d_summarized, measure = "Area", exclude = "blank", log = TRUE)

usethis::use_data(data_normalized)
