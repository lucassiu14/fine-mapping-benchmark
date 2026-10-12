# =============================================================================
# scripts/analysis/lsr_palette.R
#
# The single source of colour and shape for every LSR results figure, sourced
# by iter004_lsr_figures.R and iter004_native_cs_figure.R. A method is drawn
# with the same colour and the same shape in every figure because both
# scripts read them from here; defining them per figure is how green came to
# mean polyfun_ldsc in one figure and sbayesrc in the next.
#
# The seven are the methods displayed in the results section, in the order
# they are listed there. Okabe-Ito hues, the project's existing scheme; every
# other method is grey. Shape is a secondary encoding so identity never rests
# on colour alone - beatrice open circle, fb_xregion filled circle, so the
# original method and its extension read as a pair.
# =============================================================================
M7 <- c("susie", "polyfun_ldsc", "finemap_inf", "polyfun_oracle",
        "sbayesrc", "beatrice", "fb_xregion")
PAL7 <- c(susie          = "#0072B2",
          polyfun_ldsc   = "#009E73",
          finemap_inf    = "#56B4E9",
          polyfun_oracle = "#000000",
          sbayesrc       = "#CC79A7",
          beatrice       = "#E69F00",
          fb_xregion     = "#D55E00")
SHP7 <- c(susie = 15, polyfun_ldsc = 17, finemap_inf = 18, polyfun_oracle = 3,
          sbayesrc = 4, beatrice = 1, fb_xregion = 16)
GREY <- "#9AA0A6"            # every method outside the seven
REF  <- "#C1272D"            # reference and bound lines (y = x, y = 1 - t)
stopifnot(identical(names(PAL7), M7), identical(names(SHP7), M7),
          !anyDuplicated(PAL7), !anyDuplicated(SHP7))

# --- PDF device -------------------------------------------------------------
# The figures carry Greek and superscript characters in their axis labels (rho,
# h^2), so they were written with cairo_pdf. Homebrew R on macOS compiles the
# cairo capability in - capabilities("cairo") is TRUE - while failing to load
# libcairo at run time: ggsave then warns "failed to load cairo DLL", reports
# that it wrote the file, and leaves the old one untouched. A figure script that
# appears to succeed while changing nothing is the worst of both, and it is how
# a stale figure reaches a document.
#
# Pick a device that actually works, preferring cairo where it does. quartz is
# macOS-only but renders the same characters; the plain pdf device is the last
# resort and will mangle them, so it warns.
PDF_DEVICE <- local({
  works <- function(dev) {
    f <- tempfile(fileext = ".pdf")
    ok <- tryCatch({
      suppressWarnings(dev(filename = f, width = 2, height = 2))
      plot.new(); grDevices::dev.off()
      file.exists(f) && file.size(f) > 0
    }, error = function(e) FALSE)
    while (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(f); isTRUE(ok)
  }
  quartz_pdf <- function(filename, width, height, ...)
    grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
  if (works(grDevices::cairo_pdf)) return(grDevices::cairo_pdf)
  if (capabilities("aqua") && works(quartz_pdf)) {
    message("lsr_palette: cairo_pdf unavailable, using quartz(type = \"pdf\")")
    return(quartz_pdf)
  }
  warning("lsr_palette: neither cairo_pdf nor quartz is available; the plain pdf ",
          "device cannot render the Greek and superscript axis labels",
          call. = FALSE)
  function(filename, width, height, ...) grDevices::pdf(file = filename,
                                                        width = width, height = height)
})
