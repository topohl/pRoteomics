# House style for reader-facing review workbooks (behaviour/physiology workbook style): Arial 10, dark-blue header,
# light row rules, p values as "<0.001", values below 0.05 highlighted, a blue-white-yellow scale on signed effects,
# frozen panes, no gridlines, and a "Read me" sheet first. Needs openxlsx.

REVIEW_PALETTE <- yaml::read_yaml("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Exp9_manuscript/config/manuscript_palette.yml")
DV <- REVIEW_PALETTE$diverging

HEAD <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, fontColour = "#FFFFFF", fgFill = "#1F3D52", textDecoration = "bold",
                              halign = "center", valign = "center", wrapText = TRUE, border = "bottom", borderColour = "#FFFFFF")
LEFT <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, halign = "left", border = "bottom", borderColour = "#E5E5E5")
CTR <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, halign = "center", border = "bottom", borderColour = "#E5E5E5")
N2 <- openxlsx::createStyle(numFmt = "0.00", halign = "center")
N3 <- openxlsx::createStyle(numFmt = '[<0.001]"<0.001";0.000', halign = "center")
INT <- openxlsx::createStyle(numFmt = "0", halign = "center")
SIGP <- openxlsx::createStyle(fontColour = "#000000", bgFill = "#F2CA4E", textDecoration = "bold")
TITLE <- openxlsx::createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold", fontColour = "#1F3D52")
SUB <- openxlsx::createStyle(fontName = "Arial", fontSize = 11, textDecoration = "bold", fontColour = "#1F3D52")
WRAP <- openxlsx::createStyle(wrapText = TRUE, valign = "top")

put_table <- function(wb, sheet, df, d_cols = integer(), p_cols = integer(), sig_cols = integer(), n2_cols = integer(),
                      int_cols = integer(), left_cols = 1, freeze_col = 2, widths = NULL, lim = 1, row0 = 1) {
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df, startRow = row0, headerStyle = HEAD, withFilter = TRUE)
  nr <- nrow(df); openxlsx::setRowHeights(wb, sheet, rows = row0, heights = 34)
  widths <- if (is.null(widths)) c(28, rep(13, ncol(df) - 1)) else c(widths, rep(13, max(0, ncol(df) - length(widths))))[seq_len(ncol(df))]
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = widths)
  openxlsx::freezePane(wb, sheet, firstActiveRow = row0 + 1, firstActiveCol = freeze_col); openxlsx::showGridLines(wb, sheet, FALSE)
  if (!nr) return(invisible())
  rows <- (row0 + 1):(row0 + nr)
  openxlsx::addStyle(wb, sheet, CTR, rows = rows, cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, sheet, LEFT, rows = rows, cols = left_cols, gridExpand = TRUE, stack = TRUE)
  if (length(c(d_cols, n2_cols))) openxlsx::addStyle(wb, sheet, N2, rows = rows, cols = c(d_cols, n2_cols), gridExpand = TRUE, stack = TRUE)
  if (length(int_cols)) openxlsx::addStyle(wb, sheet, INT, rows = rows, cols = int_cols, gridExpand = TRUE, stack = TRUE)
  for (cc in d_cols) openxlsx::conditionalFormatting(wb, sheet, cols = cc, rows = rows, type = "colourScale",
                                                     style = c(DV$low, DV$mid, DV$high), rule = c(-lim, 0, lim))
  if (length(p_cols)) openxlsx::addStyle(wb, sheet, N3, rows = rows, cols = p_cols, gridExpand = TRUE, stack = TRUE)
  for (cc in sig_cols) openxlsx::conditionalFormatting(wb, sheet, cols = cc, rows = rows, type = "expression", style = SIGP,
                                                       rule = paste0("AND(ISNUMBER(", openxlsx::int2col(cc), rows[1], "),", openxlsx::int2col(cc), rows[1], "<0.05)"))
}

add_readme <- function(wb, title, blocks) {
  openxlsx::addWorksheet(wb, "Read me"); r <- 1
  openxlsx::writeData(wb, "Read me", title, startRow = r); openxlsx::addStyle(wb, "Read me", TITLE, rows = r, cols = 1); r <- r + 2
  body <- openxlsx::createStyle(fontName = "Arial", fontSize = 11, wrapText = TRUE, valign = "top")
  for (b in blocks) {
    openxlsx::writeData(wb, "Read me", b$heading, startRow = r); openxlsx::addStyle(wb, "Read me", SUB, rows = r, cols = 1); r <- r + 1
    openxlsx::writeData(wb, "Read me", data.frame(x = b$lines), startRow = r, colNames = FALSE)
    openxlsx::addStyle(wb, "Read me", body, rows = r:(r + length(b$lines) - 1), cols = 1, gridExpand = TRUE)
    r <- r + length(b$lines) + 1
  }
  openxlsx::setColWidths(wb, "Read me", cols = 1, widths = 150); openxlsx::showGridLines(wb, "Read me", FALSE)
}

col_idx <- function(df, names) match(names, names(df))
