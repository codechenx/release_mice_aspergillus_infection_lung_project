library("targets")


Sys.setenv(TAR_PROJECT = "load_and_preprocessing_data")
tar_make()

Sys.setenv(TAR_PROJECT = "fig1")
tar_make()

Sys.setenv(TAR_PROJECT = "fig2")
tar_make()

Sys.setenv(TAR_PROJECT = "fig3")
tar_make()

Sys.setenv(TAR_PROJECT = "fig4")
tar_make()

Sys.setenv(TAR_PROJECT = "fig5")
tar_make()

Sys.setenv(TAR_PROJECT = "supplementary") 
tar_make()