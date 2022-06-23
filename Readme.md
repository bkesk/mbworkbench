# Many-Body Workbench (mbworkbench)

## For Users

Instructions for using Maby-Body workbench will go here.

## For Developers

### Code Structure / Convetions

The core of many-body workbench is the `Workflow` class.
A `Workflow` instance will contain several `Block` instances where each
`Block` is a step in the workflow.
A `Workflow` instance contains the data that may be relevant throughout the workflow.
This is a design decision which allows `Block` instances to be invoked without needing
to directly communicate with one another.

Each `Block` type is defined in a subdirectory within the `/mbworkbench` directory.
Any implementation of a `Block` should be independent of other blocks and should only depend on
the presence of settings in the input yaml file, and/or on data which can be found in the
`Workflow` data.
Whenever possible, a `Block` type should implement the `check_data()` function to ensure the presence
of all necessary settings/data before the workflow attempts to run the block.

In some cases, block implementations need additional utlities beyond what is feasible to place in a single source file.
Any generally useful functions/classes should be placed within `/mbworkbench/lib/`, perhaps in a submodule.
Any code that is very specific (i.e. interfacing to/from a third-party resource) should be included
within the block implementation's directory. i.e. inside of `/mbworkbench/[block directory]/utils/` or similar.
