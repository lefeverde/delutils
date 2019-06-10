#TODO Latest RStudio makes this addin obsolete
insertAssAddin <- function() {
  # get console editor id
  # context <- rstudioapi::getConsoleEditorContext()
  # id <- context$id
  # rstudioapi::insertText("<- ", id = id)
  rstudioapi::insertText("<- ")
}

# 'foo bar'
# 'foo bar '
#
# tmp_sel <-  setSelectionRanges(document_range(c(9,9, 9,10)), id=id)
# id <- tf$id
# setCursorPosition(tmp_rng)
# tf <-  rstudioapi::getActiveDocumentContext()
#
# rstudioapi::setCursorPosition(c(9,10,9,9))
# setSelectionRanges(tf$selection[[1]]$range)
# tf$selection[[1]]
# document_range(tf$selection[[1]])
