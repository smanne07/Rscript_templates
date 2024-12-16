##New analysis project using workflowr
setwd("C:/Users/sasik/Dropbox (WherryLab)")
library(workflowr)

#Run only once
wflow_git_config(user.name = "smanne07", user.email = "sasikanthmanne@gmail.com")


wflow_start("GSE147507_tenover_RNAseq",existing = T)
wflow_build()
wflow_publish(c("analysis/index.Rmd", "analysis/about.Rmd", "analysis/license.Rmd"),
              "Publish the initial files for myproject")

wflow_git_push()

