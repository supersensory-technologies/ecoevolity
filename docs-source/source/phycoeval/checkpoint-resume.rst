|phyco_logo_long|

.. _phycocheckpoint:

##########################
Stopping and Resuming Runs
##########################

.. contents::
    :local:
    :depth: 2

Long |phyco|_ analyses can take days.
If a run is interrupted (e.g., you press ``Ctrl-C``, the machine is rebooted for
an update, or the job is killed), you do not have to start over.
When you turn on checkpointing, |phyco|_ periodically saves everything it needs
to continue the MCMC chain to a *checkpoint file*, and you can resume the run
from the checkpoint later.

**A resumed run is not just similar to an uninterrupted run; it is the same
run.**
The state, tree, and operator logs of a run that was stopped and resumed
(any number of times) are byte-for-byte identical to the logs of the same
analysis (same config file, data, and seed) that was never interrupted.

Checkpointing is off by default.
If you do not use the options below, |phyco|_ behaves exactly as before, and
its output is identical to that of previous versions.


*****************************
Turning on checkpointing
*****************************

Add ``--checkpoint`` (or ``--checkpoint-every``) to your usual command::

    phycoeval --seed 12345 --checkpoint-every 100000 phycoeval-config.yml

``--checkpoint``
    Turn on checkpointing with the default interval, which is 100 samples'
    worth of generations (i.e., 100 times ``sample_frequency`` from the
    ``mcmc_settings`` of your config file).

``--checkpoint-every N``
    Turn on checkpointing, and write a checkpoint every ``N`` generations.
    ``N`` has to be a multiple of the ``sample_frequency`` in your config file.
    A checkpoint is written right after the sample of that generation is
    logged.
    With ``--checkpoint-every 0``, no periodic checkpoints are written;
    a checkpoint is only written when you stop the run (see below).

The checkpoint is written to a file next to the log files, named like the
log files but with ``checkpoint`` in place of ``state``, ``trees``, and
``operator``; for example, the checkpoint of a run of
``phycoeval-config.yml`` is ``phycoeval-config-checkpoint-run-1.txt``
(and ``--prefix`` is added to the name as it is for the log files).
The file is replaced atomically, so there is always a complete, valid
checkpoint on disk, even if the machine loses power while a new one is being
written.
Before a checkpoint is written, all the log files are flushed to disk.
When the run finishes, the checkpoint file is deleted.

.. note::

    If you ask for checkpoints with ``--checkpoint-every``, a checkpoint is
    all you need to continue; you do **not** have to stop the run cleanly.
    If the run is killed (even with ``kill -9``) or the machine crashes, the
    most recent checkpoint is still valid.


*****************
Stopping a run
*****************

If checkpointing is on, |phyco|_ stops cleanly when it receives ``SIGINT``
(``Ctrl-C``), ``SIGTERM`` (the default signal of ``kill``, and what a system
sends when it shuts down or reboots), or ``SIGHUP``.
A signal that was already being ignored when |phyco|_ started stays ignored:
for example, ``SIGHUP`` under ``nohup``, or ``SIGINT`` for a run that a
script started in the background with ``&``.
It then

1.  finishes the generation it is in,
2.  writes a checkpoint,
3.  flushes all the log files, and
4.  prints a message that tells you how to resume the run and exits with
    status **75**.

If you send a second signal while it is doing this, |phyco|_ exits
immediately, with status **76**, without writing anything new.
The previous checkpoint is still valid.

If checkpointing is not turned on, signals are not intercepted, and they stop
|phyco|_ the way they always have.

The exit statuses of |phyco|_ are:

====== ===============================================================
Status Meaning
====== ===============================================================
0      The run finished (or ``--dry-run``).
75     The run was stopped by a signal after writing a checkpoint, and
       can be resumed with ``--resume``.
76     A second signal forced |phyco|_ to exit immediately; the last
       checkpoint written is still valid.
====== ===============================================================

Any other non-zero status (or termination by a signal) means an error.


*****************
Resuming a run
*****************

To continue a run, go to the directory you started it from, repeat the
command you used to start it, and add ``--resume`` with the path to the
checkpoint file::

    phycoeval --seed 12345 --checkpoint-every 100000 \
        --resume phycoeval-config-checkpoint-run-1.txt \
        phycoeval-config.yml

(When |phyco|_ stops, it prints this command for you.)
|phyco|_ then

*   checks that the checkpoint was written by the same analysis: the config
    file and data file have to have the same contents, and the seed, the
    ``--prefix``, the chain length, the number of threads, and the other
    options have to be the same.
    If anything is different, |phyco|_ stops with an error that lists the
    differences, and does not touch your files;
*   restores the state of the chain (the tree, the node heights and which
    nodes share them, all the other parameters, the random number generator,
    and the tuning and bookkeeping of all of the operators),
    which the checkpoint holds at full precision;
*   cuts off anything that was written to the log files after the checkpoint
    (for example, if the run was killed in the middle of writing a tree), and
    appends to the existing state, tree, and operator logs, instead of
    creating new ones;
*   carries on from the generation of the checkpoint.

When the run finishes, the tree log ends with ``END;`` like any other.
You do not have to give ``--seed`` when resuming (it is read from the
checkpoint), but if you do, it has to match.
Checkpointing stays on while resuming, and updates the checkpoint file you
resumed from.
You can add ``--dry-run`` to check that a checkpoint matches your config file
without resuming.

.. note::

    To get identical output, resume with the same build of |phyco|_
    (and so the same compiler and C++ standard library) that wrote the
    checkpoint; a checkpoint written by a build with a different C++ standard
    library is refused.
    The checkpoint format is not designed to move between versions.
    Do not edit or move the log files, and do not change the config file,
    between stopping and resuming.
    The checkpoint stores the paths of the log files as they were given
    (relative to the directory the run was started from), so resume from the
    same directory; otherwise |phyco|_ stops with an error that a log file
    does not exist, and does not touch any files.
