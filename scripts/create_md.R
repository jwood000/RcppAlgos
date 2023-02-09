dir.create("temp_md")
setwd("temp_md")

for (f in list.files("../scripts", full.names = TRUE)) {
    source(f)
    file.rename(list.files(pattern = "*_reprex.md"),
                to = paste0(gsub(".*/(.*).R", "\\1", f), ".md"))
    sapply(list.files(pattern = ".*.R"), file.remove)
    print(gsub(".*/(.*).R", "\\1", f))
}

system("mv *.md ../")
setwd("../")
system("rm -r temp_md")
system("perl -p -i -e 's/\(http(.*?)\)/\(<http$1>\)/g' *.md")
