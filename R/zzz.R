.onLoad <- function(libname, pkgname) {
  desc <- read.dcf(file.path(system.file(package = pkgname, lib.loc = libname), "DESCRIPTION"))
  packageStartupMessage("Package: ", desc[, "Package"])
  packageStartupMessage("Description: ", desc[, "Description"])
  packageStartupMessage("Version: ", desc[, "Version"])
  packageStartupMessage("Release Date: ", desc[, "Date"])
  packageStartupMessage("Authors: ", desc[, "Author"])
  packageStartupMessage("Maintainer: ", desc[, "Maintainer"])
}
