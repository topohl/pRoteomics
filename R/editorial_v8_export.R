# Editorial-v8 vector export.
#
# THE DEFECT THIS FIXES
# ---------------------
# nv_pdf() prefers rsvg::rsvg_pdf, which is vector-preserving, but falls back to
# ImageMagick at 300 dpi when rsvg is absent. rsvg is NOT installed here, so
# every final PDF in Parts 16-23 is a single full-page RGB raster: one
# /Subtype /Image of exactly 2161 x 2008 px and zero font objects. The SVG
# artwork was always fine; only the PDF conversion was wrong.
#
# THE FIX
# -------
# No SVG -> PDF conversion at all. The page is composed directly on R's
# cairo_pdf device by drawing each panel's ggplot grob into an exact millimetre
# viewport. cairo_pdf emits real text-showing operators and subset-embeds the
# font (/FontFile2, /BaseFont ...+ArialMT), so the result is genuine vector.
#
# To do that the build needs the ggplot OBJECTS, not the written SVGs. The
# frozen renderers hand their plot to nv_save_panel() and return nothing, so
# e8_capture_on() temporarily shadows nv_save_panel in the global environment:
# renderers resolve it by name at call time, so they hand their plot to the
# shadow, which records it and then delegates to the original. No frozen
# renderer is edited.

# ------------------------------------------------------------- plot capture

e8_registry <- new.env(parent = emptyenv())

e8_key <- function(id, w, h) sprintf("%s@%gx%g", id, w, h)

e8_capture_on <- function() {
  if (exists(".e8_orig_save", envir = globalenv())) return(invisible(FALSE))
  orig <- get("nv_save_panel", envir = globalenv())
  assign(".e8_orig_save", orig, envir = globalenv())
  assign("nv_save_panel", function(plot, svg_path, w_mm, h_mm) {
    assign(basename(svg_path), list(plot = plot, w = w_mm, h = h_mm),
           envir = e8_registry)
    orig(plot, svg_path, w_mm, h_mm)
  }, envir = globalenv())
  invisible(TRUE)
}

e8_capture_off <- function() {
  if (!exists(".e8_orig_save", envir = globalenv())) return(invisible(FALSE))
  assign("nv_save_panel", get(".e8_orig_save", envir = globalenv()),
         envir = globalenv())
  rm(".e8_orig_save", envir = globalenv())
  invisible(TRUE)
}

e8_registry_clear <- function() rm(list = ls(e8_registry), envir = e8_registry)

# A panel may be a plain ggplot or a patchwork composition. ggplotGrob() on a
# patchwork silently returns only part of it - that is how the Figure-2
# hippocampus drawing went missing on the first vector compose, while its
# adjacent key still rendered. Dispatch on class instead.
e8_as_grob <- function(p) {
  if (inherits(p, "patchwork")) return(patchwork::patchworkGrob(p))
  if (inherits(p, "ggplot")) return(ggplot2::ggplotGrob(p))
  if (inherits(p, "grob")) return(p)
  stop("e8_as_grob: unsupported panel class: ",
       paste(class(p), collapse = ", "), call. = FALSE)
}

# ------------------------------------------------------ vector PDF assembly

# Compose one figure page directly on cairo_pdf. svg_names maps each layout item
# to the basename the renderer wrote, so a panel used at two boxes resolves to
# the right capture.
e8_compose_pdf <- function(figure, svg_names, pdf_path) {
  W <- as.numeric(figure$width_mm)
  H <- as.numeric(figure$height_mm)
  lab_pt <- nv_pt("panel_label_pt")
  fam <- nv_palette()$typography$family

  dir_create(dirname(pdf_path))
  # cairo floors the MediaBox to whole points, which would make a 183 mm page
  # 182.7 mm. Rounding UP to the next whole point lands within 0.03 mm instead
  # of 0.28 mm short. Panel content is still placed in exact millimetres.
  pw <- ceiling(W / 25.4 * 72) / 72
  ph <- ceiling(H / 25.4 * 72) / 72
  grDevices::cairo_pdf(pdf_path, width = pw, height = ph,
                       onefile = FALSE, bg = "white", family = fam)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  for (i in seq_along(figure$layout)) {
    it <- figure$layout[[i]]
    id <- as.character(it$panel)
    x <- as.numeric(it$x); y <- as.numeric(it$y)
    w <- as.numeric(it$w); h <- as.numeric(it$h)
    nm <- svg_names[[e8_key(id, w, h)]]
    rec <- if (!is.null(nm) && exists(nm, envir = e8_registry))
      get(nm, envir = e8_registry) else NULL
    if (!is.null(rec)) {
      # grid y runs from the bottom; the contract states y from the top
      vp <- grid::viewport(x = grid::unit(x + w / 2, "mm"),
                           y = grid::unit(H - (y + h / 2), "mm"),
                           width = grid::unit(w, "mm"),
                           height = grid::unit(h, "mm"))
      grid::pushViewport(vp)
      grid::grid.draw(e8_as_grob(rec$plot))
      grid::popViewport()
    }
    lab <- as.character(it$label %||% letters[i])
    grid::grid.text(lab,
                    x = grid::unit(x, "mm"),
                    y = grid::unit(H - y - lab_pt * 0.3528 * 0.92, "mm"),
                    just = c("left", "bottom"),
                    gp = grid::gpar(fontfamily = fam, fontface = "bold",
                                    fontsize = lab_pt, col = "black"))
  }
  invisible(pdf_path)
}

# ------------------------------------------------------------- vector audit

# qpdf is the only PDF tool available here, so decompress with --qdf and read
# the object structure directly. Compressed streams cannot be grepped.
e8_pdf_probe <- function(pdf_path) {
  out <- list(page_w_mm = NA_real_, page_h_mm = NA_real_, font_count = 0L,
              embedded_fonts = 0L, n_images = 0L, largest_raster = "none",
              page_raster = FALSE, text_ops = 0L)
  if (!file.exists(pdf_path)) return(out)
  tmp <- tempfile(fileext = ".pdf")
  ok <- suppressWarnings(system2("qpdf",
    c("--qdf", "--object-streams=disable", shQuote(pdf_path), shQuote(tmp)),
    stdout = FALSE, stderr = FALSE))
  src <- if (file.exists(tmp)) tmp else pdf_path
  raw <- readBin(src, "raw", file.info(src)$size)
  s <- rawToChar(raw[raw != as.raw(0)])
  Encoding(s) <- "bytes"
  g <- function(p) regmatches(s, gregexpr(p, s, useBytes = TRUE))[[1]]

  mb <- g("/MediaBox *\\[[^]]*\\]")
  if (length(mb)) {
    v <- as.numeric(strsplit(gsub("[^0-9. ]", " ", mb[1]), " +")[[1]])
    v <- v[is.finite(v)]
    if (length(v) >= 4) {
      out$page_w_mm <- round((v[3] - v[1]) * 25.4 / 72, 1)
      out$page_h_mm <- round((v[4] - v[2]) * 25.4 / 72, 1)
    }
  }
  out$font_count <- length(g("/BaseFont"))
  out$embedded_fonts <- length(g("/FontFile[0-9]?"))
  out$text_ops <- length(g("(Tj|TJ)"))

  wid <- suppressWarnings(as.numeric(sub("/Width +", "", g("/Width +[0-9]+"))))
  hei <- suppressWarnings(as.numeric(sub("/Height +", "", g("/Height +[0-9]+"))))
  out$n_images <- length(g("/Subtype *exp/Image"))
  out$n_images <- length(g("/Subtype */Image"))
  if (length(wid) && length(hei)) {
    k <- min(length(wid), length(hei))
    area <- wid[1:k] * hei[1:k]
    j <- which.max(area)
    out$largest_raster <- sprintf("%d x %d px", wid[j], hei[j])
    # a page-sized raster is one whose pixel box matches the page at >= 150 dpi
    if (!is.na(out$page_w_mm)) {
      dpi_w <- wid[j] / (out$page_w_mm / 25.4)
      dpi_h <- hei[j] / (out$page_h_mm / 25.4)
      out$page_raster <- dpi_w >= 150 && dpi_h >= 150
    }
  }
  if (file.exists(tmp)) unlink(tmp)
  out
}

e8_vector_audit <- function(pdfs) {
  rows <- lapply(pdfs, function(p) {
    z <- e8_pdf_probe(p)
    data.frame(
      figure = sub("[.]pdf$", "", basename(p)),
      page_size_mm = sprintf("%g x %g", z$page_w_mm, z$page_h_mm),
      font_count = z$font_count,
      embedded_fonts = z$embedded_fonts,
      page_sized_raster_present = z$page_raster,
      largest_raster_dimensions = z$largest_raster,
      n_raster_objects = z$n_images,
      text_show_operators = z$text_ops,
      vector_export_pass = (!z$page_raster) && z$embedded_fonts > 0 &&
        z$text_ops > 0,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# The SAME grid composition, written to SVG via svglite. Using one routine for
# both devices guarantees the delivered SVG and PDF are the same drawing, and
# lets the composed page be previewed without a PDF rasteriser (no ghostscript
# or pdftools is available in this environment).
e8_compose_svg <- function(figure, svg_names, svg_path) {
  W <- as.numeric(figure$width_mm)
  H <- as.numeric(figure$height_mm)
  lab_pt <- nv_pt("panel_label_pt")
  fam <- nv_palette()$typography$family

  dir_create(dirname(svg_path))
  svglite::svglite(svg_path, width = W / 25.4, height = H / 25.4, bg = "white",
                   fix_text_size = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))
  for (i in seq_along(figure$layout)) {
    it <- figure$layout[[i]]
    id <- as.character(it$panel)
    x <- as.numeric(it$x); y <- as.numeric(it$y)
    w <- as.numeric(it$w); h <- as.numeric(it$h)
    nm <- svg_names[[e8_key(id, w, h)]]
    rec <- if (!is.null(nm) && exists(nm, envir = e8_registry))
      get(nm, envir = e8_registry) else NULL
    if (!is.null(rec)) {
      vp <- grid::viewport(x = grid::unit(x + w / 2, "mm"),
                           y = grid::unit(H - (y + h / 2), "mm"),
                           width = grid::unit(w, "mm"),
                           height = grid::unit(h, "mm"))
      grid::pushViewport(vp)
      grid::grid.draw(e8_as_grob(rec$plot))
      grid::popViewport()
    }
    grid::grid.text(as.character(it$label %||% letters[i]),
                    x = grid::unit(x, "mm"),
                    y = grid::unit(H - y - lab_pt * 0.3528 * 0.92, "mm"),
                    just = c("left", "bottom"),
                    gp = grid::gpar(fontfamily = fam, fontface = "bold",
                                    fontsize = lab_pt, col = "black"))
  }
  invisible(svg_path)
}
