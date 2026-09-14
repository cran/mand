# mand 3.0

## Compatibility

* Preserved the public function names and formal arguments used in
  mand 2.0 and in the accompanying book vignettes.
* Added regression tests for the public API and book vignette usage.
* Restored documentation for the 12 public functions and retained
  documentation for the package data sets.

## Dependency changes

* Removed the dependency on 'imager'.
* Reimplemented sizechange() using base R interpolation while preserving
  its published function interface.
* Moved 'caret' and 'mnormt' to optional dependencies and added
  informative checks when their functionality is requested.
* Added 'msma' as a suggested package for examples and vignettes.

## Visualization

* Replaced the previous heat-style color helper with an independent
  internal palette implemented using grDevices::colorRampPalette().
* Kept the palette helper internal and outside the public API.

## Documentation and testing

* Migrated public function documentation and NAMESPACE generation to
  roxygen2.
* Added conditional execution of examples requiring 'msma'.
* Added automated tests for public API compatibility, optional
  dependencies, image resizing, and internal color palette behavior.
* Verified successful rebuilding of all seven package vignettes.

## CRAN preparation

* Added explicit imports for functions from 'graphics', 'grDevices',
  and 'stats'.
* Declared a minimum R version of 3.5.0 for version 3 serialized data.
* Updated package metadata for the mand 3.0 book-compatible release.
