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

#include "phycoeval.hpp"


void write_phy_splash(std::ostream& out) {
    std::string v = "Version ";
    v += PROJECT_DETAILED_VERSION;
    out << string_util::banner('=') << "\n" 
        << string_util::center("Phycoeval") << "\n"
        << string_util::center("Estimating phylogenetic coevality") << "\n\n"
        << string_util::center("Part of:") << "\n"
        << string_util::center(PROJECT_NAME) << "\n"
        << string_util::center(v) << "\n"
        << string_util::banner('=') << "\n";
}

void update_log_paths(
        std::string & tree_log_path,
        std::string & state_log_path,
        std::string & operator_log_path,
        unsigned int max_number_of_attempts) {
    unsigned int ntries = 0;
    while (true) {
        if (
                (! path::exists(tree_log_path))
                && (! path::exists(state_log_path))
                && (! path::exists(operator_log_path))
            ) {
            return;
        }
        increment_log_path(tree_log_path);
        increment_log_path(state_log_path);
        increment_log_path(operator_log_path);
        ++ntries;
        if (ntries > max_number_of_attempts) {
            throw EcoevolityError("Could not generate unique output files");
        }
    }
}

void increment_log_path(
        std::string & log_path) {
    std::vector<std::string> path_elements;
    std::pair<std::string, std::string> prefix_ext;
    std::string new_suffix;
    int run_number;

    prefix_ext = path::splitext(log_path);
    path_elements = string_util::split(prefix_ext.first, '-');
    run_number = std::stoi(path_elements.back());
    path_elements.pop_back();
    ++run_number;
    
    new_suffix = "-" + std::to_string(run_number) + prefix_ext.second;

    log_path = string_util::join(path_elements, "-") + new_suffix;
}

std::string get_checkpoint_path(const std::string & state_log_path) {
    const std::string marker = "-state-run-";
    std::pair<std::string, std::string> prefix_ext = path::splitext(state_log_path);
    std::string::size_type pos = prefix_ext.first.rfind(marker);
    if (pos == std::string::npos) {
        return state_log_path + ".checkpoint.txt";
    }
    return (prefix_ext.first.substr(0, pos) +
            "-checkpoint-run-" +
            prefix_ext.first.substr(pos + marker.size()) +
            ".txt");
}

std::string get_standard_library_name() {
    std::string name = "unknown";
#if defined(_LIBCPP_VERSION)
    name = "libc++";
#elif defined(__GLIBCXX__)
    name = "libstdc++";
#endif
    return name + "/" + std::to_string(sizeof(void *) * 8) + "-bit";
}

void verify_checkpoint_identity(
        const checkpoint::CheckpointFile & file,
        const std::vector< std::pair<std::string, std::string> > & identity,
        const std::string & checkpoint_path) {
    std::map<std::string, std::string> descriptions;
    descriptions["config_hash"] = "the config file (its contents)";
    descriptions["data_hash"] = "the data file (its contents)";
    descriptions["seed"] = "the seed";
    descriptions["prefix"] = "the --prefix";
    descriptions["chain_length"] = "the chain length";
    descriptions["sample_frequency"] = "the sample frequency";
    descriptions["moves_per_generation"] = "the number of moves per generation";
    descriptions["threads"] = "the number of threads";
    descriptions["ignore_data"] = "--ignore-data";
    descriptions["strict_on_constant_sites"] = "--relax-constant-sites";
    descriptions["strict_on_missing_sites"] = "--relax-missing-sites";
    descriptions["strict_on_triallelic_sites"] = "--relax-triallelic-sites";
    descriptions["standard_library"] = "the C++ standard library the program was built with";

    std::ostringstream problems;
    for (auto const & item : identity) {
        std::string what = item.first;
        if (descriptions.count(item.first) > 0) {
            what = descriptions.at(item.first);
        }
        if (! file.has_field(item.first)) {
            problems << "    " << what << ": not recorded in the checkpoint\n";
            continue;
        }
        const std::string & saved = file.get_field(item.first);
        if (saved == item.second) {
            continue;
        }
        if (item.first == "program_version") {
            // Not necessarily a problem, but exact reproducibility is only
            // guaranteed with the same build
            std::cerr << "WARNING: The checkpoint was written by version \'"
                      << saved << "\', but this is version \'" << item.second
                      << "\'. The continued run is only guaranteed to match an "
                      << "uninterrupted run if the same build of the program "
                      << "is used.\n";
            continue;
        }
        problems << "    " << what << ": checkpoint has \'" << saved
                 << "\', but this run has \'" << item.second << "\'\n";
    }
    if (! problems.str().empty()) {
        std::ostringstream message;
        message << "\n"
                << "#######################################################################\n"
                << "###############################  ERROR  ###############################\n"
                << "Cannot resume from the checkpoint \'" << checkpoint_path << "\':\n"
                << "it was not written by a run with the same data and settings.\n"
                << "Differences:\n"
                << problems.str()
                << "To continue a run, use exactly the same config file, data, and\n"
                << "options as the run that wrote the checkpoint.\n"
                << "#######################################################################\n";
        throw EcoevolityError(message.str());
    }
}

std::string get_resume_command(
        int argc,
        char * argv[],
        const std::string & checkpoint_path) {
    auto quote = [](const std::string & s) -> std::string {
        if ((! s.empty()) &&
                (s.find_first_not_of(
                    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_@%+=:,./-")
                 == std::string::npos)) {
            return s;
        }
        std::string quoted = "\'";
        for (char c : s) {
            if (c == '\'') {
                quoted += "\'\\\'\'";
            }
            else {
                quoted += c;
            }
        }
        quoted += "\'";
        return quoted;
    };
    bool already_resuming = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "--resume") || (arg.compare(0, 9, "--resume=") == 0)) {
            already_resuming = true;
        }
    }
    std::string command = quote(argv[0]);
    if (! already_resuming) {
        command += " " + quote("--resume=" + checkpoint_path);
    }
    for (int i = 1; i < argc; ++i) {
        command += " " + quote(argv[i]);
    }
    return command;
}
