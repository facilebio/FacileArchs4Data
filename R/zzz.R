.pkgcache <- new.env(parent = emptyenv())

.onLoad <- function(libname, pgkname) {
  pkg.opts <- list(
    archs4.data_dir = NULL,
    archs4.use_cache = TRUE)

  opts <- options()
  toset <- !(names(pkg.opts) %in% names(opts))
  if (any(toset)) {
    options(pkg.opts[toset])
  }

  invisible()
}
