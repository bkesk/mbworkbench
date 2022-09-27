# Many-Body Workbench (mbworkbench)

## For Users

Instructions for using Maby-Body workbench.

### Installing


#### As a Python package

mbworkbench may be installed directly on a system as a python package using:

```bash
$ git clone https://github.com/bkesk/mbworkbench
$ pip install -r Requirements.txt
$ pip install -e .
```

mbworkbench can then be run as a module:

```bash
$ python -m mbworkbench
```

#### As a container

Alternatively, mbworkbench can be installed as a Docker / Podman container using the following syntax.
The syntax is essentially identical for docker vs podman.
Podman offers the advantage of allowing containers to be run in "rootless" mode (i.e. without requiring root
privileges at run time)

```
$ wget https://raw.githubusercontent.com/bkesk/mbworkbench/main/Containerfile
$ podman build -t mbwb:latest .
```

Rebuilding the container image may require the following to reflect recent changes

```
$ podman build --no-cache -t mbwb:latest .
```

To run the container, use the following:

```bash
podman run -ti --rm --volume $(pwd)/:/var/data/:Z mbwb:latest
```

which will start an interactive container instance where mbworkbench is
installed as a package.

Here, the `--volume $(pwd)/:/var/data/:Z` argument (or `--volume $(pwd)/:/var/data/` for Docker)
mounts the current directory as volume in the container. This allows the container to access input files 
in the current directory and allows access to output files outside of the container.
This is important since container instances are generally ephemeral.

### Basic usage

Many-body workbench is an automation tool focused specifically on "workflows"
defined as a set of "blocks" which perform a particular step in the workflow.
A workflow is specified in an input file (`input.yaml` by default) as an ordered list.

For example,

```yaml
workflow:
  - molecule
  - scf
  - write_afqmclab
```

where each item corresponds to a type of block.
In this case, the workflow consists of three steps:

1. `molecule` : a molecule is defined
2. `scf` : an scf field calculation is performed
3. `write_afqmclab` : inputs for AFQMCLab will be generated

Each of these steps may have a corresponding input block inside of 'input.yaml'.
For example, an entire input file may look like:

```yaml
geom:
  comment: "a Carbon dimer"
  atoms: |
    C 0.0 0.0  0.75
    C 0.0 0.0 -0.75
basis:
  C:
    pyscf_lib: True
    data: ccpvdz
molecule:
    spin: 0
    charge: 0
scf:
    type: rohf
write_afqmclab:
    basis: scf
    integrals: molecule
    write_wfns:
      - scf
    run_params:
      dt: 0.02
workflow:
  - molecule
  - scf
  - write_afqmclab
```

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
