#include "catch.hpp"

#include <cmath>
#include <csignal>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <streambuf>

#include "ecoevolity/checkpoint.hpp"
#include "ecoevolity/phycoeval.hpp"
#include "ecoevolity/rng.hpp"
#include "ecoevolity/path.hpp"

#include "utils_for_testing.hpp"


namespace {

std::string slurp(const std::string & path) {
    std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void spit(const std::string & path, const std::string & text) {
    std::ofstream out(path.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    out << text;
}

// Swallows everything written to it, and calls a function whenever a
// complete line of output has been written. Used to catch the moment phycoeval
// reaches some generation, by replacing the buffer of std::cout.
class LineWatcher : public std::streambuf {
    protected:
        std::string line_;
        std::function<void(const std::string &)> on_line_;

        int_type overflow(int_type c) {
            if (c == traits_type::eof()) {
                return traits_type::not_eof(c);
            }
            if (c == '\n') {
                std::string line = this->line_;
                this->line_.clear();
                this->on_line_(line);
            }
            else {
                this->line_ += static_cast<char>(c);
            }
            return c;
        }

    public:
        LineWatcher(std::function<void(const std::string &)> on_line) :
            on_line_(on_line) { }
};

class CoutRedirect {
    protected:
        std::streambuf * original_;
    public:
        CoutRedirect(std::streambuf * buffer) :
            original_(std::cout.rdbuf(buffer)) { }
        ~CoutRedirect() {
            std::cout.rdbuf(this->original_);
        }
};

int run_phycoeval(const std::vector<std::string> & arguments) {
    std::vector<std::string> storage = arguments;
    std::vector<char *> argv;
    for (auto & s : storage) {
        argv.push_back(&s[0]);
    }
    argv.push_back(NULL);
    return phycoeval_main<BasePopulationTree>(
            static_cast<int>(storage.size()),
            &argv[0]);
}

// With shared_sizes, the model is set up like many empirical analyses: one
// population size shared by all branches, a comb starting tree, a fixed
// alpha of the beta prior on node heights, and fixed state frequencies.
// Shared and fixed parameters are saved and restored differently from
// branch-specific and estimated ones.
void write_config(const std::string & path, bool shared_sizes = false) {
    std::ofstream os;
    os.open(path);
    os << "data:\n";
    os << "    ploidy: 2\n";
    os << "    constant_sites_removed: false\n";
    os << "    yaml_allele_counts:\n";
    os << "        path: species-6-genomes-2-chars-300-simulated.yml\n";
    os << "tree_model:\n";
    os << "    tree_space: generalized\n";
    os << "    starting_tree: " << (shared_sizes ? "comb" : "random") << "\n";
    os << "    tree_prior:\n";
    os << "        uniform_root_and_betas:\n";
    os << "            parameters:\n";
    os << "                root_height:\n";
    os << "                    estimate: true\n";
    os << "                    prior:\n";
    os << "                        gamma_distribution:\n";
    os << "                            shape: 10.0\n";
    os << "                            mean: 0.05\n";
    if (! shared_sizes) {
        os << "                alpha_of_node_height_beta_prior:\n";
        os << "                    value: 1.0\n";
        os << "                    estimate: true\n";
        os << "                    prior:\n";
        os << "                        exponential_distribution:\n";
        os << "                            mean: 1.0\n";
    }
    os << "branch_parameters:\n";
    os << "    population_size:\n";
    os << "        equal_population_sizes: " << (shared_sizes ? "true" : "false") << "\n";
    os << "        value: 0.004\n";
    os << "        estimate: true\n";
    os << "        prior:\n";
    os << "            gamma_distribution:\n";
    os << "                shape: 20.0\n";
    os << "                mean: 0.004\n";
    os << "mutation_parameters:\n";
    os << "    freq_1:\n";
    os << "        value: 0.5\n";
    if (shared_sizes) {
        os << "        estimate: false\n";
    }
    else {
        os << "        estimate: true\n";
        os << "        prior:\n";
        os << "            beta_distribution:\n";
        os << "                alpha: 2.0\n";
        os << "                beta: 2.0\n";
    }
    os << "    mutation_rate:\n";
    os << "        value: 1.0\n";
    os << "        estimate: false\n";
    os << "mcmc_settings:\n";
    os << "    chain_length: 120\n";
    os << "    sample_frequency: 10\n";
    os.close();
}

} // namespace


TEST_CASE("Testing checkpoint number encoding", "[checkpoint]") {

    SECTION("Testing doubles round trip exactly") {
        std::vector<double> values = {
            0.0, -0.0, 1.0, 0.1, 1.0 / 3.0, 123456789.123456789,
            std::numeric_limits<double>::min(),
            std::numeric_limits<double>::denorm_min(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
        };
        for (double x : values) {
            double y = checkpoint::string_to_double(checkpoint::double_to_string(x));
            REQUIRE(y == x);
            REQUIRE(std::signbit(y) == std::signbit(x));
        }
        REQUIRE(std::isnan(checkpoint::string_to_double(
                checkpoint::double_to_string(std::numeric_limits<double>::quiet_NaN()))));
    }

    SECTION("Testing bad numbers throw") {
        REQUIRE_THROWS_AS(checkpoint::string_to_double("1.0x"), EcoevolityError);
        REQUIRE_THROWS_AS(checkpoint::string_to_double(""), EcoevolityError);
        REQUIRE_THROWS_AS(checkpoint::string_to_unsigned("-1"), EcoevolityError);
        REQUIRE_THROWS_AS(checkpoint::string_to_unsigned("12a"), EcoevolityError);
        REQUIRE(checkpoint::string_to_integer("-12") == -12);
        REQUIRE(checkpoint::string_to_unsigned("12") == 12);
    }
}


TEST_CASE("Testing checkpoint file", "[checkpoint]") {

    SECTION("Testing round trip") {
        RandomNumberGenerator rng(123);
        std::string path = "data/tmp-checkpoint-file-" + rng.random_string(10) + ".txt";

        checkpoint::CheckpointFile file;
        file.add_field("generation", 1234);
        file.add_field("path", "some dir/with spaces.log");
        file.add_section("alpha", "line 1\nline 2\n");
        file.add_section("empty", "");
        file.write(path);
        REQUIRE(path::exists(path));
        REQUIRE(! path::exists(path + ".tmp"));

        checkpoint::CheckpointFile read = checkpoint::CheckpointFile::read(path);
        REQUIRE(read.get_unsigned("generation") == 1234);
        REQUIRE(read.get_field("path") == "some dir/with spaces.log");
        REQUIRE(read.get_section("alpha").size() == 2);
        REQUIRE(read.get_section("alpha").at(1) == "line 2");
        REQUIRE(read.get_section("empty").size() == 0);
        REQUIRE(read.to_string() == file.to_string());
        REQUIRE(! read.has_field("nope"));
        REQUIRE_THROWS_AS(read.get_field("nope"), EcoevolityError);

        // Replacing is atomic and complete
        checkpoint::CheckpointFile second;
        second.add_field("generation", 5);
        second.write(path);
        REQUIRE(checkpoint::CheckpointFile::read(path).get_unsigned("generation") == 5);

        std::remove(path.c_str());
    }

    SECTION("Testing corrupt and truncated files are detected") {
        RandomNumberGenerator rng(123);
        std::string path = "data/tmp-checkpoint-file-" + rng.random_string(10) + ".txt";
        checkpoint::CheckpointFile file;
        file.add_field("generation", 1234);
        file.add_section("alpha", "line 1\nline 2\n");
        std::string text = file.to_string();

        std::string changed = text;
        changed.replace(changed.find("1234"), 4, "1235");
        spit(path, changed);
        REQUIRE_THROWS_AS(checkpoint::CheckpointFile::read(path), EcoevolityError);

        spit(path, text.substr(0, text.size() - 10));
        REQUIRE_THROWS_AS(checkpoint::CheckpointFile::read(path), EcoevolityError);

        spit(path, "");
        REQUIRE_THROWS_AS(checkpoint::CheckpointFile::read(path), EcoevolityError);

        spit(path, text);
        REQUIRE_NOTHROW(checkpoint::CheckpointFile::read(path));
        std::remove(path.c_str());
    }
}


TEST_CASE("Testing rolling back log files", "[checkpoint]") {

    SECTION("Testing log is cut to checkpoint, and mismatches are caught") {
        RandomNumberGenerator rng(123);
        std::string path = "data/tmp-checkpoint-log-" + rng.random_string(10) + ".log";
        spit(path, "line 1\nline 2\n");
        checkpoint::CheckpointFile file;
        checkpoint::record_log_file(file, "state", path);

        spit(path, "line 1\nline 2\nline 3\nhalf a li");
        checkpoint::roll_back_log_file(file, "state", path);
        REQUIRE(slurp(path) == "line 1\nline 2\n");

        // Shorter than at the checkpoint
        spit(path, "line 1\n");
        REQUIRE_THROWS_AS(
                checkpoint::roll_back_log_file(file, "state", path),
                EcoevolityError);

        // A different file
        spit(path, "line X\nline Y\nline 3\n");
        REQUIRE_THROWS_AS(
                checkpoint::roll_back_log_file(file, "state", path),
                EcoevolityError);
        std::remove(path.c_str());
    }
}


TEST_CASE("Testing random number generator checkpointing", "[checkpoint]") {

    SECTION("Testing restored generator continues the same sequence") {
        RandomNumberGenerator rng(98765);
        // Odd number of normal draws, so that an implementation that caches
        // the second variate of a pair has one pending
        for (unsigned int i = 0; i < 101; ++i) {
            rng.uniform_real();
            rng.normal();
            rng.gamma(2.5, 0.3);
            rng.uniform_int(0, 1000);
        }
        rng.normal();
        std::ostringstream state;
        rng.write_state(state);
        std::vector<std::string> lines;
        std::istringstream in(state.str());
        std::string line;
        while (std::getline(in, line)) {
            lines.push_back(line);
        }

        RandomNumberGenerator restored(1);
        restored.read_state(lines);
        REQUIRE(restored.get_seed() == 98765);

        for (unsigned int i = 0; i < 500; ++i) {
            REQUIRE(rng.normal() == restored.normal());
            REQUIRE(rng.gamma(0.7, 1.9) == restored.gamma(0.7, 1.9));
            REQUIRE(rng.uniform_real() == restored.uniform_real());
            REQUIRE(rng.uniform_int(-5, 5000) == restored.uniform_int(-5, 5000));
            REQUIRE(rng.uniform_positive_int(77) == restored.uniform_positive_int(77));
        }
        std::vector<double> alpha = {0.5, 1.5, 3.0};
        REQUIRE(rng.dirichlet(alpha) == restored.dirichlet(alpha));
    }

    SECTION("Testing malformed state is rejected") {
        RandomNumberGenerator rng(1);
        std::vector<std::string> lines = {"seed 1", "engine 1 2 3"};
        REQUIRE_THROWS_AS(rng.read_state(lines), EcoevolityError);
    }
}


TEST_CASE("Testing phycoeval checkpoint option errors", "[checkpoint]") {

    SECTION("Testing bad options") {
        RandomNumberGenerator rng(839798);
        std::string tag = rng.random_string(10);
        std::string cfg_path = "data/tmp-config-" + tag + "-checkpoint.yml";
        write_config(cfg_path);

        // Not a multiple of the sample frequency
        REQUIRE_THROWS_AS(
                (run_phycoeval({"phycoeval", "--seed", "1", "--checkpoint-every=15", cfg_path})),
                EcoevolityError);
        REQUIRE_THROWS_AS(
                (run_phycoeval({"phycoeval", "--seed", "1", "--checkpoint-every=-10", cfg_path})),
                EcoevolityError);
        REQUIRE_THROWS_AS(
                (run_phycoeval({"phycoeval", "--seed", "1",
                        "--resume=data/tmp-no-such-checkpoint-" + tag + ".txt", cfg_path})),
                EcoevolityError);
        // Nothing was started
        REQUIRE(! path::exists("data/tmp-config-" + tag + "-checkpoint-state-run-1.log"));

        std::remove(cfg_path.c_str());
    }
}


TEST_CASE("Testing phycoeval stop and resume gives identical output",
        "[checkpoint]") {

    // Catch runs the test case once for each of these sections, each time
    // followed by the stop-and-resume test below them
    bool shared_sizes = false;
    SECTION("Testing with branch-specific population sizes") {
        shared_sizes = false;
    }
    SECTION("Testing with one shared population size") {
        shared_sizes = true;
    }

    // SIGTERM, SIGKILL-like, and resume
    {
        RandomNumberGenerator rng(839798);
        std::string tag = rng.random_string(10) + (shared_sizes ? "-shared" : "");
        const std::string cfg_path = "data/tmp-config-" + tag + "-ckpt.yml";
        const std::string base = "data/tmp-config-" + tag + "-ckpt";
        const std::string state_path = base + "-state-run-1.log";
        const std::string tree_path = base + "-trees-run-1.nex";
        const std::string oper_path = base + "-operator-run-1.log";
        const std::string checkpoint_path = base + "-checkpoint-run-1.txt";
        const std::string prefix_a = "data/tmp-a-" + tag + "-";
        write_config(cfg_path, shared_sizes);

        // Leftovers of an earlier failed run would change the output paths
        const std::string a_base = prefix_a + "tmp-config-" + tag + "-ckpt";
        const std::vector<std::string> leftovers = {
            state_path, tree_path, oper_path, checkpoint_path,
            checkpoint_path + ".tmp",
            a_base + "-state-run-1.log", a_base + "-trees-run-1.nex",
            a_base + "-operator-run-1.log"};
        for (auto const & leftover : leftovers) {
            std::remove(leftover.c_str());
        }

        const std::string seed = "20240601";

        // A: uninterrupted run, no checkpointing
        {
            LineWatcher silent([](const std::string &) { });
            CoutRedirect quiet(&silent);
            int ret = run_phycoeval({"phycoeval", "--seed", seed,
                    "--prefix", prefix_a, cfg_path});
            REQUIRE(ret == 0);
        }
        const std::string a_state = slurp(prefix_a + "tmp-config-" + tag + "-ckpt-state-run-1.log");
        const std::string a_trees = slurp(prefix_a + "tmp-config-" + tag + "-ckpt-trees-run-1.nex");
        const std::string a_oper = slurp(prefix_a + "tmp-config-" + tag + "-ckpt-operator-run-1.log");
        REQUIRE(a_trees.substr(a_trees.size() - 4) == "END;");
        REQUIRE(! path::exists(checkpoint_path));

        // B: checkpoints every 40 generations; at generation 100 (not a
        // checkpoint generation) the run is sent SIGTERM, and stops
        std::string snapshot_checkpoint;
        std::string snapshot_state;
        std::string snapshot_trees;
        std::string snapshot_oper;
        std::string snapshot_log_checkpoint;
        {
            LineWatcher watcher(
                [&](const std::string & line) {
                    if (line == "Checkpoint written after generation 80") {
                        snapshot_checkpoint = slurp(checkpoint_path);
                    }
                    if (line.compare(0, 4, "100\t") == 0) {
                        // As if the process were killed here: the checkpoint
                        // is the one from generation 80, but the logs go on
                        // to generation 100
                        snapshot_state = slurp(state_path);
                        snapshot_trees = slurp(tree_path);
                        snapshot_oper = slurp(oper_path);
                        std::raise(SIGTERM);
                    }
                });
            CoutRedirect quiet(&watcher);
            int ret = run_phycoeval({"phycoeval", "--seed", seed,
                    "--checkpoint-every=40", cfg_path});
            REQUIRE(ret == checkpoint::EXIT_STATUS_STOPPED);
        }
        REQUIRE(path::exists(checkpoint_path));
        REQUIRE(! path::exists(checkpoint_path + ".tmp"));
        {
            std::string trees = slurp(tree_path);
            REQUIRE(trees.substr(trees.size() - 4) != "END;");
        }
        REQUIRE(! snapshot_checkpoint.empty());
        // The stopped run has everything up to generation 100 and no more
        REQUIRE(slurp(state_path) == snapshot_state);

        // The stopped run is not a prefix problem: what is logged so far is
        // what the uninterrupted run logged
        REQUIRE(a_state.compare(0, snapshot_state.size(), snapshot_state) == 0);
        REQUIRE(a_trees.compare(0, snapshot_trees.size(), snapshot_trees) == 0);

        // A mismatching config is refused, and leaves the logs alone
        {
            std::string changed_cfg_path = "data/tmp-config-" + tag + "-changed.yml";
            std::string text = slurp(cfg_path);
            text.replace(text.find("chain_length: 120"), 17, "chain_length: 130");
            spit(changed_cfg_path, text);
            REQUIRE_THROWS_AS(
                    (run_phycoeval({"phycoeval", "--resume=" + checkpoint_path,
                            changed_cfg_path})),
                    EcoevolityError);
            REQUIRE(slurp(state_path) == snapshot_state);
            std::remove(changed_cfg_path.c_str());
        }

        // Resume after the SIGTERM
        {
            LineWatcher silent([](const std::string &) { });
            CoutRedirect quiet(&silent);
            int ret = run_phycoeval({"phycoeval", "--resume=" + checkpoint_path, cfg_path});
            REQUIRE(ret == 0);
        }
        REQUIRE(slurp(state_path) == a_state);
        REQUIRE(slurp(tree_path) == a_trees);
        REQUIRE(slurp(oper_path) == a_oper);
        REQUIRE(! path::exists(checkpoint_path));

        // Resume as if the run had been killed after the checkpoint at
        // generation 80, with logs that go further than the checkpoint, and
        // a half-written tree on the end
        spit(checkpoint_path, snapshot_checkpoint);
        spit(state_path, snapshot_state + "100");
        spit(tree_path, snapshot_trees + "    TREE gen1");
        spit(oper_path, snapshot_oper);
        {
            LineWatcher silent([](const std::string &) { });
            CoutRedirect quiet(&silent);
            int ret = run_phycoeval({"phycoeval", "--seed", seed,
                    "--resume=" + checkpoint_path, cfg_path});
            REQUIRE(ret == 0);
        }
        REQUIRE(slurp(state_path) == a_state);
        REQUIRE(slurp(tree_path) == a_trees);
        REQUIRE(slurp(oper_path) == a_oper);

        std::remove(cfg_path.c_str());
        std::remove(state_path.c_str());
        std::remove(tree_path.c_str());
        std::remove(oper_path.c_str());
        std::remove((prefix_a + "tmp-config-" + tag + "-ckpt-state-run-1.log").c_str());
        std::remove((prefix_a + "tmp-config-" + tag + "-ckpt-trees-run-1.nex").c_str());
        std::remove((prefix_a + "tmp-config-" + tag + "-ckpt-operator-run-1.log").c_str());
    }
}
