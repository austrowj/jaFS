source("renv/activate.R")

# Stop prompting to save workspace image.
utils::assignInNamespace(
  "q", 
  function(save = "no", status = 0, runLast = TRUE) {
    .Internal(quit(save, status, runLast))
  }, 
  "base"
)
