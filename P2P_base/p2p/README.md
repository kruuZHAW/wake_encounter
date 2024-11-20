# P2P_base

The executable model is the file `P2P` which can be found in this repository. Try to re-compile it if you have trouble running it.

Note that it is recommended to compile/run the model under *Linux*. If you have the patience, you can also try to get it to run on a different platform.

General information on `P2P` can be found in the document `P2P_documentation_ZHAW.pdf` in this repository.

## Cloning from repo
In case you don't want to re-compile `P2P` and use the executable as it is on the repository, you'll have to make it executable for Linux. You'll have to run the following command every time you get a new version from the repository:
```
chmod +x P2P
```

## Compile P2P
To re-compile the executable, make sure that you have `ifort` installed (Intel Fortran Compiler). The following instructions are to compile the source on the node `srv-lab-t-250`. Compiling it locally or on a different node might be slightly different.


Log in on `srv-lab-t-250` via SSH and set your working directory to the root of this repository.

1. Initialize the Fortran compiler:
   ```
   source /opt/intel/2021.1.1/setvars.sh
   ```
2. Compile `P2P`:
   ```
   ifort main_P2P_ZHAW.f -C -save -o P2P subroutines_ZHAW.o
   ```

## Run P2P
The documentation in `P2P_documentation_ZHAW.pdf` describes how `P2P` can be run using the shell script `while_script`. This might be a bit complicated if you are not familiar with Linux scripting. As an alternative, you can use the file `run_cases.py`. It's written in Python and uses Python plotting instead of console based plotting functionality.