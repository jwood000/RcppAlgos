# Run from the package root. pkgdown scans root Markdown files independently
# of .Rbuildignore, so keep agent instructions outside that scan while building.
with_pkgdown_sources <- function(code) {
    stopifnot(file.exists("DESCRIPTION"), dir.exists("dev"))

    agents <- normalizePath("AGENTS.md", mustWork = FALSE)
    if (file.exists(agents)) {
        backup <- tempfile(".AGENTS-", tmpdir = normalizePath("dev"),
                           fileext = ".pkgdown")
        if (!file.rename(agents, backup)) {
            stop("Could not set AGENTS.md aside for the site build.")
        }
        on.exit({
            if (file.exists(agents) || !file.rename(backup, agents)) {
                warning("Could not restore AGENTS.md; it is preserved at: ", backup)
            }
        }, add = TRUE)
    }

    # Remove old pages before pkgdown regenerates the sitemap and search index.
    if (unlink(file.path("docs", c("AGENTS.html", "AGENTS.md"))) != 0L) {
        stop("Could not remove generated AGENTS pages.")
    }

    result <- withVisible(force(code))
    fix_pkgdown_output()
    if (result$visible) result$value else invisible(result$value)
}

# Work around two pkgdown 2.2.1 output bugs without changing source documents.
fix_pkgdown_output <- function() {
    pkg <- pkgdown::as_pkgdown(".")
    not_found <- file.path(pkg$dst_path, "404.html")
    if (file.exists(not_found)) {
        html <- paste(readLines(not_found, warn = FALSE), collapse = "\n")
        fixed <- sub('<script src="[^"]+/NA">(?=\\s*window\\.MathJax\\s*=)',
                     "<script>", html, perl = TRUE)
        if (!identical(fixed, html)) writeLines(fixed, not_found)
    }

    if (isFALSE(pkg$meta[["llm-docs"]])) return(invisible())
    temp <- tempfile(fileext = ".html")
    on.exit(unlink(temp), add = TRUE)
    to <- if (rmarkdown::pandoc_available("2.10.1")) {
        "gfm+definition_lists-raw_html"
    } else {
        "gfm-raw_html"
    }
    for (path in pkgdown:::get_site_paths(pkg)) {
        src <- file.path(pkg$dst_path, path)
        dst <- sub("\\.html$", ".md", src)
        if (!file.exists(dst)) next
        html <- xml2::read_html(src)
        main <- xml2::xml_find_first(html, ".//main")
        if (inherits(main, "xml_missing") ||
            !grepl("\\\\[([]|\\$\\$", xml2::xml_text(main))) next

        # Keep pkgdown's Markdown formatting, but let Pandoc recognize TeX.
        pkgdown:::simplify_page_header(main)
        pkgdown:::simplify_anchors(main)
        pkgdown:::simplify_code(main)
        pkgdown:::simplify_popovers_to_footnotes(main)
        pkgdown:::simplify_lifecycle_badges(main)
        pkgdown:::simplify_dls(main)
        pkgdown:::create_absolute_links(main, pkgdown:::full_url(pkg, path))
        xml2::write_html(main, temp, format = FALSE)
        rmarkdown::pandoc_convert(
            input = temp, output = dst,
            from = "html+tex_math_dollars+tex_math_single_backslash",
            to = to
        )
    }

    # llms.txt combines the Markdown indexes, which can also contain math.
    indexes <- file.path(pkg$dst_path,
                         c("index.md", "reference/index.md", "articles/index.md"))
    if (file.exists(file.path(pkg$dst_path, "llms.txt"))) {
        text <- lapply(indexes, pkgdown:::read_file_if_exists)
        writeLines(c(text[[1]], "", text[[2]], "", text[[3]]),
                   file.path(pkg$dst_path, "llms.txt"))
    }
    invisible()
}
