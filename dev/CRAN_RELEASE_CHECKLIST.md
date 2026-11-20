# CRAN Release Checklist (RcppAlgos)

This checklist is optimized for small patch releases and CRAN pretests
(CITATION readability, URL validity, clean tarball contents).

---

## 0) Before you start

- [ ] Confirm you are on a clean release branch (no WIP for the next large release).
- [ ] `git status` is clean (or only includes intentional release edits).
- [ ] Review staged changes (ignore pkgdown output as needed):
  - [ ] `git diff --stat`
  - [ ] `git diff <prev_tag>..HEAD -- . ':(exclude)docs'`

---

## 1) Versioning + changelog

- [ ] Bump version in `DESCRIPTION` (e.g., `2.9.4` → `2.9.5`).
- [ ] Match the version constant in `src/version.cpp` to `DESCRIPTION` and
      rebuild: package loading checks that the DLL and package versions agree.
- [ ] Update `NEWS.md` with a new top entry.
- [ ] Update `inst/NEWS.Rd` to match `NEWS.md`.
- [ ] Keep unreleased-version markers consistent between both NEWS files.

---

## 2) CITATION (CRAN incoming pretest requirement)

- [ ] Update `inst/CITATION`.
- [ ] Ensure `inst/CITATION` contains **no dynamic calls** that require the package
      to be installed/available (e.g., no `packageVersion()`, `Sys.Date()`,
      `utils::packageDescription()`, `citation()`, etc.).
- [ ] Ensure URLs are canonical (e.g., trailing slashes where applicable).
- [ ] Local verification:
  - [ ] `R -q -e "tools::checkRdaFiles('.')"` (optional)
  - [ ] Build + check tarball (see section 5) and confirm no:
    - “Reading CITATION file fails … when package is not installed.”

---

## 3) URLs + metadata

- [ ] Run URL checks (fix true errors):
  - [ ] `R -q -e "urlchecker::url_check()"`
- [ ] Fix:
  - [ ] Invalid fragments (bad `#...` URLs)
  - [ ] True 404s
  - [ ] Permanent redirects (301) by switching to canonical URL (often add trailing slash)
- [ ] Confirm `DESCRIPTION` includes canonical URLs (pkgdown + GitHub), e.g.:
  - [ ] `URL: https://jwood000.github.io/RcppAlgos/, https://github.com/jwood000/RcppAlgos`
  - [ ] `BugReports: https://github.com/jwood000/RcppAlgos/issues`

---

## 4) Build hygiene (.Rbuildignore)

- [ ] Ensure non-package files are excluded in `.Rbuildignore`, at minimum:
  - [ ] `^docs$` (pkgdown output)
  - [ ] `^_pkgdown\.yml$`
  - [ ] `^dev$` (development files, including `dev/CRAN_RELEASE_CHECKLIST.md`)
  - [ ] Any local tooling directories (e.g., `^ai_git_helpers$`)
- [ ] Confirm `R CMD build` tarball does not include website/config extras.

---

## 5) Local checks (must be clean)

Run checks from a clean state.

- [ ] `R CMD build .`
- [ ] `R CMD check RcppAlgos_<version>.tar.gz --as-cran`
- [ ] Run the optional native regressions, which are excluded from package checks:
      `Rscript --vanilla .github/development/native-regressions/run.R`.
- [ ] Review reverse-dependency results against the last published CRAN version;
      include changes inherited from any intervening unreleased version.

Verify:
- [ ] No ERRORs or WARNINGs.
- [ ] Any NOTEs are understood/acceptable (avoid new NOTES if possible).
- [ ] In particular:
  - [ ] CITATION readable in incoming feasibility checks
  - [ ] No invalid URLs (especially malformed GitHub links)
  - [ ] No new “non-standard file/directory found at top level”

---

## 6) pkgdown site (optional for CRAN, recommended for releases)

If maintaining the site for each release:

- [ ] Update `_pkgdown.yml` as needed (bootstrap, url).
- [ ] Run the site build through `dev/build_site.R` (see quick commands below).
- [ ] Confirm no generated `AGENTS` pages or search/sitemap entries remain.
- [ ] Commit site output under `docs/` (if that is the repo policy).
- [ ] If pkgdown warns about missing alt text, consider fixing in README/vignettes.

---

## 7) External builders

- [ ] Win-builder:
  - [ ] devel
  - [ ] release
  - [ ] oldrelease
  - [ ] Review overall check time, including installation and compilation,
        against any timing feedback from CRAN; test timings alone are insufficient.
- [ ] (Recommended) CI matrix / GitHub Actions:
  - [ ] at least Linux + macOS
  - [ ] include r-devel somewhere if practical

---

## 8) Final review + tagging

- [ ] Review final diff (excluding docs):
  - [ ] `git diff <prev_tag>..HEAD -- . ':(exclude)docs'`
- [ ] Confirm `DESCRIPTION` fields (at least):
  - [ ] `Version`
  - [ ] `URL` (canonical)
  - [ ] `BugReports`
- [ ] Create release tag:
  - [ ] `git tag -a v<version> -m "CRAN patch release <version>"`
  - [ ] `git push --tags`

---

## 9) CRAN submission

- [ ] Prepare `cran-comments.md` (if used), include:
  - [ ] Summary of changes
  - [ ] Any NOTE explanations (if any remain)
  - [ ] Explain how each issue from the previous submission was addressed,
        including timing changes and new platform measurements when relevant
- [ ] Submit via `devtools::submit_cran()` or CRAN web form.
- [ ] For resubmissions, put the response to CRAN's feedback in the web form's
      optional comment field and accept the confirmation email.

See the [CRAN submission and resubmission policy](https://cran.r-project.org/web/packages/policies.html#Submission)
for the tarball, R-devel check, reverse-dependency, and submission requirements.

---

## Quick commands (copy/paste)

### Diff ignoring pkgdown output

```sh
git diff <old> <new> -- . ':(exclude)docs'
```

### Build + check

```sh
R CMD build .
R CMD check RcppAlgos_*.tar.gz --as-cran
```

### URL check

```r
urlchecker::url_check()
```

### pkgdown

```r
source("dev/build_site.R")
with_pkgdown_sources(pkgdown::build_site())
```

The helper temporarily sets `AGENTS.md` aside and restores it when the build
finishes or fails. It also removes stale generated AGENTS pages before pkgdown
rebuilds the search index and sitemap. A plain `pkgdown::build_site()` includes
root Markdown files even when `.Rbuildignore` excludes them from the package.

## Docker UBSan (r-devel-ubsan-clang) sanity check (recommended for C++ changes)

This is a quick way to catch undefined behavior (UB) that may not show up in
regular checks. Uses the rocker image: `rocker/r-devel-ubsan-clang`.

### 1) Start Docker
- [ ] Open Docker Desktop and ensure it is running.

### 2) Run the container (mount repo)
> `--rm` removes the container when you exit; the repo is mounted in.

```sh
docker run --platform linux/amd64 -m=10g --cap-add SYS_PTRACE \
  --name=r-devel-ubsan-clang \
  -v ~/RcppAlgos:/RSource/RcppAlgos \
  --rm -ti rocker/r-devel-ubsan-clang /bin/bash
```

### 3) Install system deps (inside container)

```sh
apt-get update
apt-get install -y libgmp3-dev
```

### 4) Install R deps (inside container)

```sh
R -e 'install.packages(c("cpp11","gmp","testthat","microbenchmark","numbers","partitions"), repos="https://cloud.r-project.org")'
R -e 'install.packages("RcppBigIntAlgos", repos="https://cloud.r-project.org")'
```

### 5) Build + check package (inside container)

```sh
R CMD build /RSource/RcppAlgos
Rdevel CMD check --no-examples --no-vignettes RcppAlgos_*.tar.gz
```

Notes:

* Use `--no-examples --no-vignettes` for speed. Remove these flags when debugging
  vignette/example-related issues.
* If you are investigating memory/UB issues deeply, consider running without
  `--no-vignettes` to widen coverage.

---

## Debugging UBSan traces (addr2line)

Sometimes UBSan prints stack frames with offsets into the shared object, e.g.:

`/RcppAlgos.Rcheck/RcppAlgos/libs/RcppAlgos.so+0x3034db`

You can map these offsets back to source lines using `addr2line`.

### 1) Locate the compiled shared object and the offset

Example:

* SO: `/RcppAlgos.Rcheck/RcppAlgos/libs/RcppAlgos.so`
* Offset: `0x3034db`

### 2) Convert offset to file:line

```sh
addr2line -e /RcppAlgos.Rcheck/RcppAlgos/libs/RcppAlgos.so -f -C 0x3034db
```

Recommended:

* `-f` prints function name
* `-C` demangles C++ symbols

### 3) If symbolization is missing / warnings appear

If you see messages like:
`WARNING: Failed to use and restart external symbolizer!`

Try installing/ensuring LLVM symbolizer tools are present, then re-run.
(Exact package names can vary across base images.)
