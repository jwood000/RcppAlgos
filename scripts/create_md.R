dir.create("/Users/owacasa/RcppAlgos/temp_md",
           showWarnings = FALSE, recursive = TRUE)
setwd("/Users/owacasa/RcppAlgos/temp_md")

scripts <- list.files("../scripts", pattern = "\\.R$", full.names = TRUE)
scripts <- scripts[!grepl("create_md\\b", basename(scripts))]

for (f in scripts) {
    if (!dir.exists(f)) {
        source(f)
        print(f)
        file.rename(list.files(pattern = "*_reprex.md"),
                    to = paste0(gsub(".*/(.*).R", "\\1", f), ".md"))
        sapply(list.files(pattern = ".*.R"), file.remove)
        sapply(list.files(pattern = ".*.out"), file.remove)
        print(gsub(".*/(.*).R", "\\1", f))
    }
}

system("mv *.md ../")
setwd("../")
system("rm -r temp_md")
## Run on command line... the escapes are killing me!
## perl -p -i -e 's{\((https?://.*?)\)}{(<$1>)}g' *.md
##
## This command removes the additional coding blocks:
##
## code...
## ```
##
## ``` r
##
## more code
##
## -- to --
##
## code
##
## more code
##
## perl -0777p -i -e 's/```\n\n``` r\n//g' *.md
##
## Update the date! Change to the vignette dir and run the following:
## perl -pi -e 's/^date:\s*".*"/date: "'"$(date +%Y-%m-%d)"'"/' *.Rmd

