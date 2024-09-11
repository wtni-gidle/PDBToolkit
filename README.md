# PDBToolkit

`PDBToolkit` is a Python library for protein structure handling, featuring modules for PDB operations and CASP-related workflows.

## Setup

Before using `PDBToolkit`, please follow these steps:

### Configure `PYTHONPATH`

To ensure Python can locate PDBToolkit, set the `PYTHONPATH` environment variable:

1. Edit your shell configuration file (e.g., `~/.bashrc`)

    ```shell
    vim ~/.bashrc
    ```

2. Add the following line to the end of the file:

    ```shell
    export PYTHONPATH=/path/to/this/repository:$PYTHONPATH
    ```

3. Reload the shell configuration:

    ```shell
    source ~/.bashrc
    ```

### Update Configuration File

Modify the paths in the [`config.py`](PDBToolkit/config.py) file to ensure that the library can locate the necessary tools for proper functionality.

## TODO

- [ ] write jobs in batch and automaticly submit
- [ ] write af3 jobs
