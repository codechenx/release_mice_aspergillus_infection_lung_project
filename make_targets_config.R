library(targets)

tar_config_set(
  script = "0.load_and_preprocessing_data.R",
  store = "cache/0.load_and_preprocessing_data",
  project = "load_and_preprocessing_data",
  use_crew=TRUE,
  workers = 20,
  config = "_targets.yaml"
)

tar_config_set(
  script = "1.fig1.R",
  store = "cache/1.fig1",
  project = "fig1",
  use_crew=TRUE,
  workers = 20,
  config = "_targets.yaml"
)

tar_config_set(
  script = "2.fig2.R",
  store = "cache/2.fig2",
  project = "fig2",
  use_crew=TRUE,
  workers = 20,
  config = "_targets.yaml"
)

tar_config_set(
  script = "3.fig3.R",
  store = "cache/3.fig3",
  project = "fig3",
  use_crew=TRUE,
  workers = 20,
  config = "_targets.yaml"
)

tar_config_set(
  script = "4.fig4.R",
  store = "cache/4.fig4",
  project = "fig4",
  use_crew=TRUE,
  workers = 20,
  config = "_targets.yaml"
)

tar_config_set(
  script = "5.fig5.R",
  store = "cache/5.fig5",
  project = "fig5",
  use_crew=TRUE,
  workers = 20,
  config = "_targets.yaml"
)


tar_config_set(
  script = "99.supplementary.R",
  store = "cache/99.supplementary",
  project = "supplementary",
  use_crew=TRUE,
  workers = 1,
  config = "_targets.yaml"
)