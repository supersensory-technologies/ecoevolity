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

#ifndef ECOEVOLITY_MCMC_HPP
#define ECOEVOLITY_MCMC_HPP

#include "assert.hpp"
#include "error.hpp"
#include "checkpoint.hpp"

#include "tree.hpp"
#include "general_tree_operator.hpp"
#include "general_tree_operator_schedule.hpp"

/**
 * Write a checkpoint from which the run can be continued, with output
 * identical to an uninterrupted run. All log files are flushed to disk first,
 * and their sizes recorded, so that anything logged after the checkpoint can
 * be cut off when resuming.
 */
template<class TreeType>
inline void write_mcmc_checkpoint(
        const checkpoint::Control & control,
        const unsigned int generation,
        const unsigned int gen_of_last_state_log,
        const unsigned int gen_of_last_operator_log,
        const RandomNumberGenerator & rng,
        const TreeType & tree,
        const GeneralTreeOperatorSchedule<TreeType> & operator_schedule,
        std::ostream & tree_log_stream,
        std::ostream & state_log_stream,
        std::ostream & operator_log_stream) {
    tree_log_stream.flush();
    state_log_stream.flush();
    operator_log_stream.flush();

    checkpoint::CheckpointFile file;
    file.add_field("format", checkpoint::FORMAT_VERSION);
    for (auto const & item : control.identity) {
        file.add_field(item.first, item.second);
    }
    file.add_field("checkpoint_every", control.every);
    file.add_field("generation", generation);
    file.add_field("gen_of_last_state_log", gen_of_last_state_log);
    file.add_field("gen_of_last_operator_log", gen_of_last_operator_log);
    // Used to verify that the restored state is the state that was saved
    file.add_field("ln_likelihood",
            checkpoint::double_to_string(tree.get_log_likelihood_value()));
    file.add_field("ln_prior",
            checkpoint::double_to_string(tree.get_log_prior_density_value()));
    checkpoint::record_log_file(file, "tree", control.tree_log_path);
    checkpoint::record_log_file(file, "state", control.state_log_path);
    checkpoint::record_log_file(file, "operator", control.operator_log_path);

    std::ostringstream rng_state;
    rng.write_state(rng_state);
    file.add_section("rng", rng_state.str());
    std::ostringstream tree_state;
    tree.write_checkpoint_state(tree_state);
    file.add_section("tree", tree_state.str());
    std::ostringstream operator_state;
    operator_schedule.write_checkpoint_state(operator_state);
    file.add_section("operators", operator_state.str());

    file.write(control.path);
}

/**
 * Run the MCMC.
 *
 * If checkpoint_control is null or not enabled, the chain runs from start to
 * finish and writes what it always has. Otherwise, checkpoints are written at
 * the requested intervals, and the run stops cleanly (after writing a
 * checkpoint) when SIGINT or SIGTERM is received. The log streams are assumed
 * to be positioned at the end of the checkpointed output when resuming.
 *
 * @return  True if the chain ran to the end; false if it was stopped by a
 *          signal (see checkpoint_control).
 */
template<class TreeType>
inline bool mcmc(
        RandomNumberGenerator & rng,
        TreeType & tree,
        GeneralTreeOperatorSchedule<TreeType> & operator_schedule,
        const unsigned int chain_length,
        const unsigned int sample_frequency,
        const unsigned int number_of_moves_per_generation,
        std::ostream & tree_log_stream,
        std::ostream & state_log_stream,
        std::ostream & operator_log_stream,
        std::ostream & std_output_stream,
        const std::string & logging_delimiter = "\t",
        const unsigned int logging_precision = 18,
        const unsigned int nthreads = 1,
        checkpoint::Control * checkpoint_control = nullptr) {
    const bool checkpointing = (checkpoint_control && checkpoint_control->enabled);
    const bool resuming = (checkpointing && checkpoint_control->resuming);

    tree_log_stream.precision(logging_precision);
    state_log_stream.precision(logging_precision);
    operator_log_stream.precision(logging_precision);

    unsigned int start_gen = 0;
    unsigned int gen_of_last_state_log = 0;
    unsigned int gen_of_last_operator_log = 0;

    if (resuming) {
        const checkpoint::CheckpointFile & saved = checkpoint_control->resume_file;
        rng.read_state(saved.get_section("rng"));
        tree.read_checkpoint_state(saved.get_section("tree"));
        operator_schedule.read_checkpoint_state(saved.get_section("operators"));
        start_gen = saved.get_unsigned("generation");
        gen_of_last_state_log = saved.get_unsigned("gen_of_last_state_log");
        gen_of_last_operator_log = saved.get_unsigned("gen_of_last_operator_log");
        if (start_gen > chain_length) {
            throw EcoevolityError(
                    "The checkpoint is from a later generation than the chain length");
        }
        tree.write_state_log_header(std_output_stream, logging_delimiter, true);
    }
    else {
        tree.write_state_log_header(state_log_stream, logging_delimiter);
        tree.write_state_log_header(std_output_stream, logging_delimiter, true);

        tree_log_stream << "#NEXUS" << std::endl;
        tree.write_nexus_taxa_block(tree_log_stream);
        tree_log_stream << "\nBEGIN TREES;\n";
    }

    tree.make_dirty();
    tree.compute_log_likelihood_and_prior(nthreads);
    if (tree.get_log_likelihood_value() == -std::numeric_limits<double>::infinity()) {
        std::ostringstream message;
        message << "\n\n"
                << "#######################################################################\n"
                << "###############################  ERROR  ###############################\n"
                << "The initial model state has a probability density of zero.\n"
                << "#######################################################################\n";
        throw EcoevolityError(message.str());
    }
    if (std::isnan(tree.get_log_likelihood_value())) {
        std::ostringstream message;
        message << "\n\n"
                << "#######################################################################\n"
                << "###############################  ERROR  ###############################\n"
                << "The initial model state has a NAN log likelihood.\n"
                << "#######################################################################\n";
        throw EcoevolityError(message.str());
    }
    if (resuming) {
        // The likelihood and prior are functions of the restored state, so
        // they have to be exactly what they were when the checkpoint was
        // written. If not, the state was not restored faithfully, and the
        // run would silently diverge from the one that was stopped.
        const checkpoint::CheckpointFile & saved = checkpoint_control->resume_file;
        if ((tree.get_log_likelihood_value() !=
                checkpoint::string_to_double(saved.get_field("ln_likelihood"))) ||
            (tree.get_log_prior_density_value() !=
                checkpoint::string_to_double(saved.get_field("ln_prior")))) {
            throw EcoevolityError(
                    "The log likelihood or log prior of the state restored from the "
                    "checkpoint differs from the values saved in it, so the checkpoint "
                    "cannot be resumed exactly");
        }
        std_output_stream << "Resuming MCMC after generation " << start_gen << "\n";
        tree.log_state(std_output_stream, start_gen, logging_delimiter, true);
    }
    else {
        tree.log_state(state_log_stream, 0, logging_delimiter);
        tree.log_state(std_output_stream, 0, logging_delimiter, true);
        tree.log_nexus_tree(tree_log_stream, 0, true, logging_precision);
    }

    std::shared_ptr< GeneralTreeOperatorTemplate< TreeType > > op;
    unsigned int gen;
    for (gen = start_gen; gen < chain_length; ++gen) {
        for (unsigned int move_count = 0;
                move_count < number_of_moves_per_generation;
                ++move_count) {
            op = operator_schedule.draw_operator(rng);
            op->operate_with_helpers(rng,
                    &tree,
                    nthreads, 1, 1);
        }

        if ((gen + 1) % sample_frequency == 0) {
            tree.log_state(state_log_stream, gen + 1, logging_delimiter);
            tree.log_nexus_tree(tree_log_stream, gen + 1, true, logging_precision);
            gen_of_last_state_log = gen;
            // Log every 10th sample to std out
            if ((gen + 1) % (sample_frequency * 10) == 0) {
                tree.log_state(std_output_stream, gen + 1, logging_delimiter, true);
                // Log operator performance every 100 samples
                if ((gen + 1) % (sample_frequency * 100) == 0) {
                    operator_log_stream << "generation " << gen + 1 << ":\n";
                    operator_schedule.write_operator_rates(operator_log_stream);
                    gen_of_last_operator_log = gen;
                }
            }
        }

        if (checkpointing) {
            bool checkpoint_is_current = false;
            if ((checkpoint_control->every > 0) &&
                    (((gen + 1) % checkpoint_control->every) == 0)) {
                write_mcmc_checkpoint<TreeType>(
                        *checkpoint_control,
                        gen + 1,
                        gen_of_last_state_log,
                        gen_of_last_operator_log,
                        rng,
                        tree,
                        operator_schedule,
                        tree_log_stream,
                        state_log_stream,
                        operator_log_stream);
                checkpoint_is_current = true;
                std_output_stream << "Checkpoint written after generation "
                        << gen + 1 << std::endl;
            }
            // Stop if asked to, unless the chain is done anyway (in which
            // case finishing only takes a moment)
            int stop_signal = checkpoint::pending_stop_signal();
            if ((stop_signal != 0) && ((gen + 1) < chain_length)) {
                if (! checkpoint_is_current) {
                    write_mcmc_checkpoint<TreeType>(
                            *checkpoint_control,
                            gen + 1,
                            gen_of_last_state_log,
                            gen_of_last_operator_log,
                            rng,
                            tree,
                            operator_schedule,
                            tree_log_stream,
                            state_log_stream,
                            operator_log_stream);
                }
                checkpoint_control->stopped_by_signal = stop_signal;
                checkpoint_control->stopped_at_generation = gen + 1;
                return false;
            }
        }
    }
    // Make sure last generation is reported
    if (gen > (gen_of_last_state_log + 1)) {
        tree.log_state(state_log_stream, gen + 1, logging_delimiter);
        tree.log_state(std_output_stream, gen + 1, logging_delimiter, true);
        tree.log_nexus_tree(tree_log_stream, gen + 1, true, logging_precision);
    }
    if (gen > (gen_of_last_operator_log + 1)) {
        operator_log_stream << "generation " << gen + 1 << ":\n";
        operator_schedule.write_operator_rates(operator_log_stream);
    }
    tree_log_stream << "END;";
    std_output_stream << "\nOperator stats:\n";
    operator_schedule.write_operator_rates(std_output_stream);
    std_output_stream << "\n";
    return true;
}

#endif
