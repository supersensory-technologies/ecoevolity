This directory contains some tools to aid in developing and testing the
ecoevolity code. Most users of ecoevolity can safely ignore this directory.

# Running the ecoevolity test suite on a cluster

I was having issues building the test suite on the compute nodes of our cluster
(bus errors despite the build working under an identical environment on the
head node).

I created the scripts in this directory to make testing more cluster-friendly.
To use them:

1.  Run the `prep-tests.sh` script on the head node. It will simply build the
    executable for the ecoevolity test suite; this runs quickly without using
    much comp resources.

    If it is run without any arguments, the script will build the ecoevolity
    tests in `../test-suite-build`. If arguments are given, the first arg is
    used as the directory where the tests are built.

    This script will also write the location of the test-suite executable to a
    file `./env-test-location.sh`. The next script will use this info.

    Before building the tests, the script will source some files if they exist:

    -   `modules-to-load.sh`: This file is meant to hold cluster-specific
        module loading commands. But, I suppose it can contain any environment
        setup.

    -   `../dependencies/env-*.sh`: The `../dependencies` dir has a few scripts
        for downloading and building some of ecoevolity's dependencies, and
        storing their locations in `env-*.sh` files. Most of these dependencies
        are bundled with ecoevolity, but building them with the scripts in the
        `dependencies` directory will allow the `prep-tests.sh` script to build
        the test suite much faster, because it doesn't have to build the
        bundled dependency code. So, if you are going to be running the tests
        regularly, it's nice to have the dependencies pre-built to save some
        time.

2.  Submit the `\.run-tests.sh` script to your job scheduler to run on a
    compute node. E.g.,

        sbatch -p -N 1 -n 1 -t 150:00:00 run-tests.sh

    The `run-tests.sh` will run the ecoevolity test suite executable built by
    `prep-tests.sh`. It gets the location of the executable from the
    `env-test-location.sh` file written by the `prept-tests.sh` script.

    Before running the tests it will source all of the same files as
    `prep-tests.sh` to make sure any libraries dynamically linked to test-suite
    executable can be found.

# Testing phycoeval checkpointing and resuming

The `test-checkpoint-resume.sh` script is an end-to-end test of phycoeval's
`--checkpoint-every`, `--checkpoint`, and `--resume` options. It uses the
synthetic data and config in `../test/data` (`config-checkpoint-resume.yml`).
It checks that, for a fixed seed, the state, tree, and operator logs of

-   an uninterrupted run,
-   a run stopped with `SIGTERM` (and again with `SIGINT`) and resumed, and
-   a run killed with `SIGKILL` after a periodic checkpoint and resumed

are byte-identical, and (optionally) that a run without the new options is
byte-identical to the output of an unmodified build. To do the latter, build the
upstream `v1.1.1` tag in a separate build directory, and give the script both
executables:

    dev-tools/test-checkpoint-resume.sh path/to/build/src/ecoevolity/phycoeval path/to/v1.1.1-build/src/ecoevolity/phycoeval

The script takes a few minutes, and exits with a non-zero status if any check
fails. The same behavior is also covered by `test_checkpoint.cpp` in the unit
test suite.
