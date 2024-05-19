find_program(SPHINX_EXEC NAMES sphinx-build)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Sphinx
                                  "Failed to find sphinx-build executable"
                                  SPHINX_EXEC)
