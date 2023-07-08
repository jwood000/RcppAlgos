dir.create("/Users/josephwood/RcppAlgos/temp_md")
setwd("/Users/josephwood/RcppAlgos/temp_md")

my_scripts <- list.files("../scripts", full.names = TRUE)
my_scripts <- my_scripts[which(!grepl("create_md", my_scripts))]

for (f in my_scripts) {
    source(f)
    print(f)
    file.rename(list.files(pattern = "*_reprex.md"),
                to = paste0(gsub(".*/(.*).R", "\\1", f), ".md"))
    sapply(list.files(pattern = ".*.R"), file.remove)
    sapply(list.files(pattern = ".*.out"), file.remove)
    print(gsub(".*/(.*).R", "\\1", f))
}

system("mv *.md ../")
setwd("../")
system("rm -r temp_md")
## Run on command line... the escapes are killing me!
## perl -p -i -e 's/\(http(.*?)\)/\(<http$1>\)/g' *.md
