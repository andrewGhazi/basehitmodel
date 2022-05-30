.onAttach = function(libname, pkgname) {
  ver = utils::packageVersion("basehitmodel")
  packageStartupMessage("This is basehitmodel version ", ver)
  packageStartupMessage("- Parallelize: Before calling basehitmodel::model_proteins_separately(), run future::plan() in a way that's appropriate for your system.")
  packageStartupMessage("- Show progress: Before calling basehitmodel::model_proteins_separately(), run library(progressr); handlers(global=TRUE)")
}
