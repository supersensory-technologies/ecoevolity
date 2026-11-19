/******************************************************************************
 * Copyright (C) 2015-2016 Jamie R. Oaks.
 *
 * This file is part of Ecoevolity.
 *
 * Ecoevolity is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Ecoevolity is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with Ecoevolity.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#ifndef PHYCOEVAL_HPP
#define PHYCOEVAL_HPP

#include <limits>
#include <cstdio>
#include <fstream>
#include <utility>
#include <time.h>

#include "cpp-optparse/OptionParser.h"

#include "version.hpp"
#include "error.hpp"
#include "rng.hpp"
#include "path.hpp"
#include "string_util.hpp"
#include "general_tree_settings.hpp"
#include "settings_io.hpp"
#include "mcmc.hpp"
#include "checkpoint.hpp"


void write_phy_splash(std::ostream& out);

void update_log_paths(
        std::string & tree_log_path,
        std::string & state_log_path,
        std::string & operator_log_path,
        unsigned int max_number_of_attempts = 10000);

void increment_log_path(std::string & log_path);

// Path of the checkpoint file that goes with a run, derived from the path of
// its state log (e.g., "x-state-run-2.log" gives "x-checkpoint-run-2.txt")
std::string get_checkpoint_path(const std::string & state_log_path);

// Name of the C++ standard library the program was built with, because
// random number generation (and so a resumed run) depends on it
std::string get_standard_library_name();

// Throw an informative error if a checkpoint was not written by a run with
// the given settings
void verify_checkpoint_identity(
        const checkpoint::CheckpointFile & file,
        const std::vector< std::pair<std::string, std::string> > & identity,
        const std::string & checkpoint_path);

// The command line to continue a run from a checkpoint: the command that
// was run, with the option to resume added
std::string get_resume_command(
        int argc,
        char * argv[],
        const std::string & checkpoint_path);

template <class TreeType>
int phycoeval_main(int argc, char * argv[]) {

    write_phy_splash(std::cout);
    std::cout << "\n";

    const std::string usage = 
        "usage: %prog [OPTIONS] YAML-CONFIG-FILE";

    std::ostringstream version_ss;
    version_ss << PROJECT_NAME << " version " << PROJECT_DETAILED_VERSION;
    const std::string version = version_ss.str();

    const std::string description =
        "Phycoeval: Estimating phylogenetic coevality";
    // const std::string epilog =
    //     "Epilog goes here...";

    optparse::OptionParser parser = optparse::OptionParser()
            .usage(usage)
            .version(version)
            .description(description);
            // .epilog(epilog);

    parser.set_defaults("ignore_data", "false");

    parser.add_option("--seed")
            .action("store")
            .type("long")
            .dest("seed")
            .help("Seed for random number generator. Default: Set from clock.");
    parser.add_option("--ignore-data")
            .action("store_true")
            .dest("ignore_data")
            .help("Ignore data to sample from the prior distribution. Default: "
                  "Use data to sample from the posterior distribution");
#ifdef BUILD_WITH_THREADS
    parser.add_option("--nthreads")
            .action("store")
            .type("unsigned int")
            .dest("nthreads")
            .set_default("1")
            .help("Number of threads to use for likelihood calculations. "
                  "Default: 1 (no multithreading). If you are using "
                  "the \'--ignore-data\' option, no likelihood calculations "
                  "will be performed, and so no multithreading is used.");
#endif
    parser.add_option("--prefix")
            .action("store")
            .dest("prefix")
            .set_default("")
            .help("Optional string to prefix all output files.");
    parser.add_option("--relax-constant-sites")
            .action("store_true")
            .dest("relax_constant_sites")
            .help("By default, if you specify \'constant_sites_removed = true\' "
                  "and constant sites are found, phycoeval throws an error. "
                  "With this option, phycoeval will automatically ignore the "
                  "constant sites and only issue a warning (and correct for "
                  "constant sites in the likelihood calculation). Please make sure "
                  "you understand what you are doing when you use this option."
                );
    parser.add_option("--relax-missing-sites")
            .action("store_true")
            .dest("relax_missing_sites")
            .help("By default, if a column is found for which there is no data "
                  "across all populations, phycoeval throws an error. "
                  "With this option, phycoeval will automatically ignore such "
                  "sites and only issue a warning."
                );
    parser.add_option("--relax-triallelic-sites")
            .action("store_true")
            .dest("relax_triallelic_sites")
            .help("By default, if a DNA site is found for which there is more "
                  "than two nucleotide states, phycoeval throws an error. "
                  "With this option, phycoeval will automatically recode such "
                  "sites as biallelic and only issue a warning. These sites "
                  "are recoded by assigning state 0 to the first nucleotide "
                  "found and state 1 to all others. If you do not wish to "
                  "recode such sites and prefer to ignore them, please remove "
                  "all sites with more than two nucleotide states from your "
                  "DNA alignments. NOTE: only alignments of nucleotides are "
                  "affected by this option, not alignments of standard "
                  "characters (i.e., 0, 1, 2)."
                );
    parser.add_option("--dry-run")
            .action("store_true")
            .dest("dry_run")
            .help("Do not run analysis; only process and report settings.");
    parser.add_option("--checkpoint")
            .action("store_true")
            .dest("checkpoint")
            .help("Turn on checkpointing: write a checkpoint file every "
                  "100 samples (or at the interval given by "
                  "\'--checkpoint-every\'), and, if the run is stopped with "
                  "SIGINT (Ctrl-C) or SIGTERM, finish the current generation, "
                  "write a checkpoint, and exit with status "
                  "75. A run can be continued from its checkpoint with "
                  "\'--resume\'; with the same build of phycoeval, the "
                  "output is identical to that of a run that was never "
                  "interrupted. A second signal exits "
                  "immediately (status 76). Default: No checkpointing, "
                  "signals are not handled.");
    parser.add_option("--checkpoint-every")
            .action("store")
            .type("long")
            .dest("checkpoint_every")
            .help("Write a checkpoint every N generations (turns on "
                  "checkpointing, as does \'--checkpoint\'). N must be a "
                  "multiple of the sample frequency. With N = 0, a "
                  "checkpoint is only written when the run is stopped by a "
                  "signal. Default: 100 samples worth of generations, if "
                  "checkpointing is on.");
    parser.add_option("--resume")
            .action("store")
            .dest("resume")
            .help("Continue a run from the given checkpoint file. Use the "
                  "same config file, data, and options as the run that wrote "
                  "the checkpoint (the seed, if given, must match), and run "
                  "it from the same directory, because the log paths are "
                  "stored as they were given. The "
                  "state, tree, and operator logs of the run are appended "
                  "to. Checkpointing stays on, and the checkpoint file is "
                  "updated in place.");

    optparse::Values& options = parser.parse_args(argc, argv);
    std::vector<std::string> args = parser.args();

    checkpoint::Control checkpoint_control;
    std::string resume_path;
    const bool resume_requested = options.is_set_by_user("resume");
    if (resume_requested) {
        resume_path = options.get("resume").get_str();
        if (! path::exists(resume_path)) {
            throw EcoevolityError("Checkpoint file \'" + resume_path +
                    "\' does not exist");
        }
        checkpoint_control.resume_file = checkpoint::CheckpointFile::read(resume_path);
        checkpoint_control.resuming = true;
        if (checkpoint_control.resume_file.get_unsigned("format") !=
                checkpoint::FORMAT_VERSION) {
            throw EcoevolityError("Checkpoint file \'" + resume_path +
                    "\' has an unsupported format version");
        }
    }

    RandomNumberGenerator rng;
    long seed_opt;
    if (resume_requested) {
        seed_opt = static_cast<long>(checkpoint::string_to_integer(
                checkpoint_control.resume_file.get_field("seed")));
        if (options.is_set_by_user("seed")) {
            long user_seed = options.get("seed");
            if (user_seed != seed_opt) {
                std::ostringstream message;
                message << "The seed (" << user_seed << ") does not match the "
                        << "seed of the run that wrote the checkpoint ("
                        << seed_opt << ")";
                throw EcoevolityError(message.str());
            }
        }
    }
    else if (options.is_set_by_user("seed")) {
        seed_opt = options.get("seed");
    }
    else {
        seed_opt = rng.uniform_int(1, std::numeric_limits<int>::max() - 1);
    }
    const long seed = seed_opt;
    rng.set_seed(seed);
    std::cout << "Seed: " << seed << std::endl;

    const bool dry_run = options.get("dry_run");

    const bool ignore_data = options.get("ignore_data");
    if (ignore_data) {
        std::cout << "Ignoring data in order to sample from the prior distribution..." << std::endl;
    }
    else {
        std::cout << "Using data in order to sample from the posterior distribution..." << std::endl;
    }

    const bool strict_on_constant_sites = (! options.get("relax_constant_sites"));
    const bool strict_on_missing_sites = (! options.get("relax_missing_sites"));
    const bool strict_on_triallelic_sites = (! options.get("relax_triallelic_sites"));

#ifdef BUILD_WITH_THREADS 
    unsigned int nthreads = options.get("nthreads");
#else
    unsigned int nthreads = 1;
#endif

    if (args.size() < 1) {
        throw EcoevolityError("Path to YAML-formatted config file is required");
    }
    if (args.size() > 1) {
        throw EcoevolityError("Too many arguments; only one config file is allowed");
    }
    const std::string config_path = args.at(0);
    if (! path::exists(config_path)) {
        throw EcoevolityError("Config file \'" + config_path +
                "\' does not exist");
    }
    if (! path::isfile(config_path)) {
        throw EcoevolityError("Config path \'" + config_path +
                "\' is not a regular file");
    }
    std::cout << "Config path: " << config_path << std::endl;

    std::cout << "Parsing config file..." << std::endl;
    PopulationTreeSettings settings(config_path);

    std::cout << "Configuring model..." << std::endl;

    TreeType tree(
            settings,
            rng,
            strict_on_constant_sites,
            strict_on_missing_sites,
            strict_on_triallelic_sites,
            false // store_seq_loci_info
            );

    GeneralTreeOperatorSchedule<BasePopulationTree> operator_schedule(
            settings.operator_settings, tree.get_leaf_node_count());

    std::cout << "\n" << string_util::banner('-') << "\n";
    write_settings(std::cout, settings, operator_schedule);
    std::cout << string_util::banner('-') << "\n\n";

    unsigned int n_moves_per_generation = tree.get_leaf_node_count();

    if (ignore_data) {
        tree.ignore_data();
    }
    else {
        tree.use_data();
    }

    unsigned int logging_precision = 18;

    std::string tree_log_path = settings.get_tree_log_path();
    std::string state_log_path = settings.get_state_log_path();
    std::string operator_log_path = settings.get_operator_log_path();
    if (options.is_set_by_user("prefix")) {
        std::string output_prefix = options.get("prefix").get_str();
        tree_log_path = output_prefix + path::basename(tree_log_path);
        state_log_path = output_prefix + path::basename(state_log_path);
        operator_log_path = output_prefix + path::basename(operator_log_path);
    }
    if (resume_requested) {
        // Continue appending to the logs of the run that wrote the checkpoint
        tree_log_path = checkpoint_control.resume_file.get_field("log_tree_path");
        state_log_path = checkpoint_control.resume_file.get_field("log_state_path");
        operator_log_path = checkpoint_control.resume_file.get_field("log_operator_path");
    }
    else {
        update_log_paths(tree_log_path, state_log_path, operator_log_path);
    }

    const bool checkpoint_flag = options.get("checkpoint");
    long checkpoint_every_opt = 0;
    if (options.is_set_by_user("checkpoint_every")) {
        checkpoint_every_opt = options.get("checkpoint_every");
        if (checkpoint_every_opt < 0) {
            throw EcoevolityError("\'--checkpoint-every\' cannot be negative");
        }
    }
    if (checkpoint_flag || (checkpoint_every_opt > 0) || resume_requested) {
        const unsigned long sample_freq = settings.get_sample_frequency();
        checkpoint_control.enabled = true;
        if (options.is_set_by_user("checkpoint_every")) {
            checkpoint_control.every = checkpoint_every_opt;
        }
        else if (resume_requested) {
            checkpoint_control.every = checkpoint_control.resume_file.get_unsigned(
                    "checkpoint_every");
        }
        else {
            checkpoint_control.every = sample_freq * 100;
        }
        if ((checkpoint_control.every % sample_freq) != 0) {
            std::ostringstream message;
            message << "The checkpoint interval (" << checkpoint_control.every
                    << ") must be a multiple of the sample frequency ("
                    << sample_freq << ")";
            throw EcoevolityError(message.str());
        }
        checkpoint_control.tree_log_path = tree_log_path;
        checkpoint_control.state_log_path = state_log_path;
        checkpoint_control.operator_log_path = operator_log_path;
        if (resume_requested) {
            checkpoint_control.path = resume_path;
        }
        else {
            checkpoint_control.path = get_checkpoint_path(state_log_path);
        }

        // Everything that has to be the same for a continued run to be the
        // run that was stopped
        std::string prefix_opt = "";
        if (options.is_set_by_user("prefix")) {
            prefix_opt = options.get("prefix").get_str();
        }
        std::vector< std::pair<std::string, std::string> > & identity =
                checkpoint_control.identity;
        identity.push_back(std::make_pair(std::string("program_version"),
                std::string(PROJECT_DETAILED_VERSION)));
        identity.push_back(std::make_pair(std::string("standard_library"),
                get_standard_library_name()));
        identity.push_back(std::make_pair(std::string("config_hash"),
                checkpoint::hash_file(config_path)));
        identity.push_back(std::make_pair(std::string("data_hash"),
                checkpoint::hash_file(settings.data_settings.get_path())));
        identity.push_back(std::make_pair(std::string("seed"),
                std::to_string(seed)));
        identity.push_back(std::make_pair(std::string("prefix"), prefix_opt));
        identity.push_back(std::make_pair(std::string("chain_length"),
                std::to_string(settings.get_chain_length())));
        identity.push_back(std::make_pair(std::string("sample_frequency"),
                std::to_string(settings.get_sample_frequency())));
        identity.push_back(std::make_pair(std::string("moves_per_generation"),
                std::to_string(n_moves_per_generation)));
        identity.push_back(std::make_pair(std::string("threads"),
                std::to_string(nthreads)));
        identity.push_back(std::make_pair(std::string("ignore_data"),
                std::string(ignore_data ? "1" : "0")));
        identity.push_back(std::make_pair(std::string("strict_on_constant_sites"),
                std::string(strict_on_constant_sites ? "1" : "0")));
        identity.push_back(std::make_pair(std::string("strict_on_missing_sites"),
                std::string(strict_on_missing_sites ? "1" : "0")));
        identity.push_back(std::make_pair(std::string("strict_on_triallelic_sites"),
                std::string(strict_on_triallelic_sites ? "1" : "0")));
        if (resume_requested) {
            verify_checkpoint_identity(
                    checkpoint_control.resume_file,
                    identity,
                    resume_path);
        }
    }

    std::cout << "\n" << string_util::banner('-') << "\n";
    tree.write_data_summary(std::cout);
    std::cout << string_util::banner('-') << "\n\n";

    std::cout << "Number of threads: " << nthreads << std::endl;

    if (dry_run) {
        return 0;
    }

    if (! resume_requested) {
        if (path::exists(tree_log_path)) {
            std::ostringstream message;
            message << "ERROR: The tree log file \'"
                    << tree_log_path
                    << "\' already exists!\n";
            throw EcoevolityError(message.str());
        }
        if (path::exists(state_log_path)) {
            std::ostringstream message;
            message << "ERROR: The state log file \'"
                    << state_log_path
                    << "\' already exists!\n";
            throw EcoevolityError(message.str());
        }
        if (path::exists(operator_log_path)) {
            std::ostringstream message;
            message << "ERROR: The operator log file \'"
                    << operator_log_path
                    << "\' already exists!\n";
            throw EcoevolityError(message.str());
        }
        if (checkpoint_control.enabled && path::exists(checkpoint_control.path)) {
            std::ostringstream message;
            message << "ERROR: The checkpoint file \'"
                    << checkpoint_control.path
                    << "\' already exists!\n";
            throw EcoevolityError(message.str());
        }
    }

    std::ofstream tree_log_stream;
    std::ofstream state_log_stream;
    std::ofstream operator_log_stream;

    if (resume_requested) {
        // Check that all the logs are the ones the checkpoint was taken
        // from, before changing any of them
        checkpoint::check_log_file(
                checkpoint_control.resume_file, "tree", tree_log_path);
        checkpoint::check_log_file(
                checkpoint_control.resume_file, "state", state_log_path);
        checkpoint::check_log_file(
                checkpoint_control.resume_file, "operator", operator_log_path);
        // Cut off anything logged after the checkpoint, so that we continue
        // appending exactly where the checkpointed run left off
        checkpoint::roll_back_log_file(
                checkpoint_control.resume_file, "tree", tree_log_path);
        checkpoint::roll_back_log_file(
                checkpoint_control.resume_file, "state", state_log_path);
        checkpoint::roll_back_log_file(
                checkpoint_control.resume_file, "operator", operator_log_path);
        tree_log_stream.open(tree_log_path, std::ios::out | std::ios::app);
        state_log_stream.open(state_log_path, std::ios::out | std::ios::app);
        operator_log_stream.open(operator_log_path, std::ios::out | std::ios::app);
    }
    else {
        tree_log_stream.open(tree_log_path);
        state_log_stream.open(state_log_path);
        operator_log_stream.open(operator_log_path);
    }
    
    if (! tree_log_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not open tree log file \'"
                << tree_log_path
                << "\'\n";
        throw EcoevolityError(message.str());
    }
    if (! state_log_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not open state log file \'"
                << state_log_path
                << "\'\n";
        throw EcoevolityError(message.str());
    }
    if (! operator_log_stream.is_open()) {
        std::ostringstream message;
        message << "ERROR: Could not open operator log file \'"
                << operator_log_path
                << "\'\n";
        throw EcoevolityError(message.str());
    }

    std::cout << "Tree log path: " << tree_log_path << std::endl;
    std::cout << "State log path: " << state_log_path << std::endl;
    std::cout << "Operator log path: " << operator_log_path << std::endl;

    time_t start;
    time_t finish;
    time(&start);

    if (resume_requested) {
        std::cout << "Resuming MCMC from checkpoint: " << resume_path << std::endl;
    }
    else {
        std::cout << "Firing up MCMC..." << std::endl;
    }
    if (checkpoint_control.enabled) {
        std::cout << "Checkpoint path: " << checkpoint_control.path << "\n";
        std::cout << "Checkpoint interval: ";
        if (checkpoint_control.every > 0) {
            std::cout << checkpoint_control.every << " generations" << std::endl;
        }
        else {
            std::cout << "only when stopped by a signal" << std::endl;
        }
    }
    bool completed = true;
    {
        // Signals are only intercepted if checkpointing is turned on
        checkpoint::StopSignalGuard signal_guard(checkpoint_control.enabled);
        completed = mcmc<TreeType>(
                rng,
                tree,
                operator_schedule,
                settings.get_chain_length(),
                settings.get_sample_frequency(),
                n_moves_per_generation,
                tree_log_stream,
                state_log_stream,
                operator_log_stream,
                std::cout,
                "\t",
                logging_precision,
                nthreads,
                checkpoint_control.enabled ? &checkpoint_control : nullptr);
    }

    tree_log_stream.close();
    state_log_stream.close();
    operator_log_stream.close();

    time(&finish);
    double duration = difftime(finish, start);
    std::cout << "Runtime: " << duration << " seconds." << std::endl;

    if (! completed) {
        const std::string resume_command = get_resume_command(
                argc,
                argv,
                checkpoint_control.path);
        std::ostringstream message;
        message << "\n"
                << "Stopped by signal " << checkpoint_control.stopped_by_signal
                << " after generation " << checkpoint_control.stopped_at_generation
                << " of " << settings.get_chain_length() << ".\n"
                << "Checkpoint written to: " << checkpoint_control.path << "\n"
                << "The logs are consistent with the checkpoint. To continue "
                << "the run, repeat the command with the option\n"
                << "\'--resume=" << checkpoint_control.path << "\' added:\n"
                << "    " << resume_command << "\n"
                << "Exit status: " << checkpoint::EXIT_STATUS_STOPPED
                << " (stopped; can be resumed).\n";
        std::cout << message.str() << std::flush;
        std::cerr << message.str() << std::flush;
        return checkpoint::EXIT_STATUS_STOPPED;
    }

    if (checkpoint_control.enabled) {
        // The run is complete, so the checkpoint is no longer needed
        std::remove(checkpoint_control.path.c_str());
        std::remove((checkpoint_control.path + ".tmp").c_str());
    }

    return 0;
}

#endif
