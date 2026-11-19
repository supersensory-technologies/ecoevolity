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

// Low-level utilities for checkpointing and resuming an MCMC run:
//
//   - exact (hex-float) text encoding of doubles
//   - a small, checksummed text file format for checkpoints
//   - atomic, durable file writing
//   - bookkeeping needed to roll log files back to a checkpoint
//   - SIGINT/SIGTERM handling for a graceful stop

#ifndef ECOEVOLITY_CHECKPOINT_HPP
#define ECOEVOLITY_CHECKPOINT_HPP

#include <atomic>
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "error.hpp"

namespace checkpoint {

// Version of the checkpoint file format
const unsigned int FORMAT_VERSION = 1;

// Exit status when a run stopped on a signal after writing a checkpoint, and
// can be continued with --resume
const int EXIT_STATUS_STOPPED = 75;
// Exit status when a second signal forced an immediate exit; the most recent
// checkpoint (if any) is still valid, but nothing new was written
const int EXIT_STATUS_FORCED = 76;

// Number of bytes at the end of a log file that are hashed to verify that a
// log file really is the one a checkpoint was taken from
const unsigned long TAIL_BYTES = 65536;


//////////////////////////////////////////////////////////////////////////////
// Exact text encoding of numbers
//////////////////////////////////////////////////////////////////////////////

// Hexadecimal floating point (e.g., 0x1.8p+1) round-trips every double
// exactly, including infinities and NaN
inline std::string double_to_string(double x) {
    char buffer[64];
    int n = std::snprintf(buffer, sizeof(buffer), "%a", x);
    if ((n < 0) || (n >= static_cast<int>(sizeof(buffer)))) {
        throw EcoevolityError("checkpoint: could not encode a floating point number");
    }
    return std::string(buffer);
}

inline double string_to_double(const std::string & s) {
    if (s.empty()) {
        throw EcoevolityError("checkpoint: empty floating point field");
    }
    char * end = NULL;
    double x = std::strtod(s.c_str(), &end);
    if ((end == s.c_str()) || (*end != '\0')) {
        throw EcoevolityError(
                "checkpoint: could not parse floating point field \'" + s + "\'");
    }
    return x;
}

inline unsigned long long string_to_unsigned(const std::string & s) {
    if (s.empty() || (s.find_first_not_of("0123456789") != std::string::npos)) {
        throw EcoevolityError(
                "checkpoint: could not parse unsigned integer field \'" + s + "\'");
    }
    errno = 0;
    unsigned long long x = std::strtoull(s.c_str(), NULL, 10);
    if (errno != 0) {
        throw EcoevolityError(
                "checkpoint: unsigned integer field out of range \'" + s + "\'");
    }
    return x;
}

inline long long string_to_integer(const std::string & s) {
    if (s.empty()) {
        throw EcoevolityError("checkpoint: empty integer field");
    }
    if (s.at(0) == '-') {
        return -static_cast<long long>(string_to_unsigned(s.substr(1)));
    }
    return static_cast<long long>(string_to_unsigned(s));
}

inline bool string_to_bool(const std::string & s) {
    if (s == "1") {
        return true;
    }
    if (s == "0") {
        return false;
    }
    throw EcoevolityError(
            "checkpoint: could not parse boolean field \'" + s + "\'");
}

inline std::vector<std::string> split_tokens(const std::string & line) {
    std::vector<std::string> tokens;
    std::istringstream stream(line);
    std::string token;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}


//////////////////////////////////////////////////////////////////////////////
// Hashing (FNV-1a, 64 bit); used to detect mismatched or corrupted files, not
// for security
//////////////////////////////////////////////////////////////////////////////

class Hash {
    protected:
        unsigned long long state_ = 14695981039346656037ULL;

    public:
        void update(const char * data, size_t length) {
            for (size_t i = 0; i < length; ++i) {
                this->state_ ^= static_cast<unsigned char>(data[i]);
                this->state_ *= 1099511628211ULL;
            }
        }
        void update(const std::string & s) {
            this->update(s.data(), s.size());
        }
        std::string hex() const {
            char buffer[32];
            std::snprintf(buffer, sizeof(buffer), "%016llx", this->state_);
            return std::string(buffer);
        }
};

inline std::string hash_string(const std::string & s) {
    Hash h;
    h.update(s);
    return h.hex();
}

// Hash of all bytes of a file
inline std::string hash_file(const std::string & path) {
    std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
    if (! in.is_open()) {
        throw EcoevolityError(
                "checkpoint: could not open \'" + path + "\' to hash it");
    }
    Hash h;
    std::vector<char> buffer(1 << 20);
    while (in) {
        in.read(&buffer[0], static_cast<std::streamsize>(buffer.size()));
        h.update(&buffer[0], static_cast<size_t>(in.gcount()));
    }
    if (in.bad()) {
        throw EcoevolityError("checkpoint: error reading \'" + path + "\'");
    }
    return h.hex();
}


//////////////////////////////////////////////////////////////////////////////
// Files
//////////////////////////////////////////////////////////////////////////////

inline unsigned long long file_size(const std::string & path) {
    struct stat info;
    if (::stat(path.c_str(), &info) != 0) {
        throw EcoevolityError(
                "checkpoint: could not get size of \'" + path + "\': " +
                std::strerror(errno));
    }
    return static_cast<unsigned long long>(info.st_size);
}

// Force the contents of a file (or directory) out to stable storage
inline bool sync_path(const std::string & path) {
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        return false;
    }
    int status = -1;
#ifdef F_FULLFSYNC
    status = ::fcntl(fd, F_FULLFSYNC);
#endif
    if (status != 0) {
        status = ::fsync(fd);
    }
    ::close(fd);
    return (status == 0);
}

inline void sync_parent_directory(const std::string & path) {
    std::string::size_type pos = path.find_last_of('/');
    std::string dir;
    if (pos == std::string::npos) {
        dir = ".";
    }
    else if (pos == 0) {
        dir = "/";
    }
    else {
        dir = path.substr(0, pos);
    }
    // Not every file system supports this; the file itself was already synced
    sync_path(dir);
}

// Replace the file at path with content such that a reader (or a crash)
// always sees either the complete old file or the complete new file
inline void write_file_atomically(
        const std::string & path,
        const std::string & content) {
    const std::string tmp_path = path + ".tmp";
    {
        std::ofstream out(tmp_path.c_str(),
                std::ios::out | std::ios::binary | std::ios::trunc);
        if (! out.is_open()) {
            throw EcoevolityError(
                    "checkpoint: could not open \'" + tmp_path + "\' for writing");
        }
        out.write(content.data(), static_cast<std::streamsize>(content.size()));
        out.flush();
        if (! out) {
            throw EcoevolityError(
                    "checkpoint: error writing \'" + tmp_path + "\'");
        }
    }
    if (! sync_path(tmp_path)) {
        throw EcoevolityError(
                "checkpoint: could not sync \'" + tmp_path + "\' to disk");
    }
    if (std::rename(tmp_path.c_str(), path.c_str()) != 0) {
        throw EcoevolityError(
                "checkpoint: could not rename \'" + tmp_path + "\' to \'" +
                path + "\': " + std::strerror(errno));
    }
    sync_parent_directory(path);
}

inline void truncate_file(const std::string & path, unsigned long long size) {
    if (::truncate(path.c_str(), static_cast<off_t>(size)) != 0) {
        throw EcoevolityError(
                "checkpoint: could not truncate \'" + path + "\': " +
                std::strerror(errno));
    }
}

// Hash of the last (at most) TAIL_BYTES bytes before the given offset in a
// file; tail_length is set to the number of bytes hashed
inline std::string hash_file_tail(
        const std::string & path,
        unsigned long long end_offset,
        unsigned long & tail_length) {
    std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
    if (! in.is_open()) {
        throw EcoevolityError(
                "checkpoint: could not open \'" + path + "\' to hash it");
    }
    tail_length = TAIL_BYTES;
    if (end_offset < tail_length) {
        tail_length = static_cast<unsigned long>(end_offset);
    }
    in.seekg(static_cast<std::streamoff>(end_offset - tail_length));
    std::vector<char> buffer(tail_length + 1);
    in.read(&buffer[0], static_cast<std::streamsize>(tail_length));
    if (static_cast<unsigned long>(in.gcount()) != tail_length) {
        throw EcoevolityError(
                "checkpoint: could not read the end of \'" + path + "\'");
    }
    Hash h;
    h.update(&buffer[0], tail_length);
    return h.hex();
}


//////////////////////////////////////////////////////////////////////////////
// Checkpoint file
//
// Layout (plain text, one item per line):
//
//     ecoevolity-checkpoint
//     <key>: <value>                (any number, no newline in either part)
//     @begin <section>
//     <lines>                       (must not start with '@')
//     @end <section>
//     @checksum <hash of all preceding bytes>
//////////////////////////////////////////////////////////////////////////////

class CheckpointFile {
    protected:
        std::vector< std::pair<std::string, std::string> > fields_;
        std::vector< std::pair<std::string, std::vector<std::string> > > sections_;

    public:
        static const char * magic() { return "ecoevolity-checkpoint"; }

        CheckpointFile() { }

        void add_field(const std::string & key, const std::string & value) {
            if ((key.empty()) ||
                    (key.find_first_of(":\n ") != std::string::npos) ||
                    (value.find('\n') != std::string::npos)) {
                throw EcoevolityError(
                        "checkpoint: invalid field \'" + key + "\'");
            }
            this->fields_.push_back(std::make_pair(key, value));
        }
        void add_field(const std::string & key, unsigned long long value) {
            this->add_field(key, std::to_string(value));
        }

        void add_section(
                const std::string & name,
                const std::string & text) {
            std::vector<std::string> lines;
            std::istringstream stream(text);
            std::string line;
            while (std::getline(stream, line)) {
                if ((! line.empty()) && (line.at(0) == '@')) {
                    throw EcoevolityError(
                            "checkpoint: invalid line in section \'" + name + "\'");
                }
                lines.push_back(line);
            }
            this->sections_.push_back(std::make_pair(name, lines));
        }

        bool has_field(const std::string & key) const {
            for (auto const & kv : this->fields_) {
                if (kv.first == key) {
                    return true;
                }
            }
            return false;
        }
        const std::string & get_field(const std::string & key) const {
            for (auto const & kv : this->fields_) {
                if (kv.first == key) {
                    return kv.second;
                }
            }
            throw EcoevolityError(
                    "checkpoint: file has no field \'" + key + "\'");
        }
        unsigned long long get_unsigned(const std::string & key) const {
            return string_to_unsigned(this->get_field(key));
        }

        const std::vector<std::string> & get_section(
                const std::string & name) const {
            for (auto const & kv : this->sections_) {
                if (kv.first == name) {
                    return kv.second;
                }
            }
            throw EcoevolityError(
                    "checkpoint: file has no section \'" + name + "\'");
        }

        std::string to_string() const {
            std::ostringstream out;
            out << magic() << "\n";
            for (auto const & kv : this->fields_) {
                out << kv.first << ": " << kv.second << "\n";
            }
            for (auto const & section : this->sections_) {
                out << "@begin " << section.first << "\n";
                for (auto const & line : section.second) {
                    out << line << "\n";
                }
                out << "@end " << section.first << "\n";
            }
            const std::string body = out.str();
            return body + "@checksum " + hash_string(body) + "\n";
        }

        void write(const std::string & path) const {
            write_file_atomically(path, this->to_string());
        }

        static CheckpointFile read(const std::string & path) {
            std::ifstream in(path.c_str(), std::ios::in | std::ios::binary);
            if (! in.is_open()) {
                throw EcoevolityError(
                        "Could not open checkpoint file \'" + path + "\'");
            }
            std::ostringstream buffer;
            buffer << in.rdbuf();
            const std::string text = buffer.str();
            const std::string corrupt =
                    "Checkpoint file \'" + path + "\' is corrupt or incomplete";

            // Verify the checksum, which is the last line
            if ((text.size() < 2) || (text.at(text.size() - 1) != '\n')) {
                throw EcoevolityError(corrupt);
            }
            std::string::size_type last_line_start =
                    text.find_last_of('\n', text.size() - 2);
            if (last_line_start == std::string::npos) {
                throw EcoevolityError(corrupt);
            }
            ++last_line_start;
            const std::string body = text.substr(0, last_line_start);
            const std::string trailer = text.substr(
                    last_line_start,
                    text.size() - last_line_start - 1);
            const std::string prefix = "@checksum ";
            if ((trailer.compare(0, prefix.size(), prefix) != 0) ||
                    (trailer.substr(prefix.size()) != hash_string(body))) {
                throw EcoevolityError(corrupt);
            }

            CheckpointFile file;
            std::istringstream stream(body);
            std::string line;
            if ((! std::getline(stream, line)) || (line != magic())) {
                throw EcoevolityError(
                        "File \'" + path + "\' is not a checkpoint file");
            }
            std::string section_name;
            std::string section_text;
            bool in_section = false;
            while (std::getline(stream, line)) {
                if (in_section) {
                    if (line == "@end " + section_name) {
                        file.add_section(section_name, section_text);
                        in_section = false;
                    }
                    else {
                        section_text += line + "\n";
                    }
                    continue;
                }
                if (line.compare(0, 7, "@begin ") == 0) {
                    section_name = line.substr(7);
                    section_text.clear();
                    in_section = true;
                    continue;
                }
                std::string::size_type pos = line.find(": ");
                if (pos == std::string::npos) {
                    throw EcoevolityError(corrupt);
                }
                file.add_field(line.substr(0, pos), line.substr(pos + 2));
            }
            if (in_section) {
                throw EcoevolityError(corrupt);
            }
            return file;
        }
};


//////////////////////////////////////////////////////////////////////////////
// Log files
//
// A checkpoint records the size of each log file (and a hash of its tail) at
// the moment the checkpoint was taken. To resume, anything written after
// that is cut off, so that the continued run appends exactly what the
// original run would have written next.
//////////////////////////////////////////////////////////////////////////////

// Flush a log file to disk and record where it ends
inline void record_log_file(
        CheckpointFile & file,
        const std::string & tag,
        const std::string & path) {
    if (! sync_path(path)) {
        throw EcoevolityError(
                "checkpoint: could not sync \'" + path + "\' to disk");
    }
    unsigned long long size = file_size(path);
    unsigned long tail_length = 0;
    std::string tail_hash = hash_file_tail(path, size, tail_length);
    file.add_field("log_" + tag + "_path", path);
    std::ostringstream value;
    value << size << " " << tail_length << " " << tail_hash;
    file.add_field("log_" + tag, value.str());
}

// Make sure a log file is the one the checkpoint was taken from; returns the
// size it had when the checkpoint was written
inline unsigned long long check_log_file(
        const CheckpointFile & file,
        const std::string & tag,
        const std::string & path) {
    std::vector<std::string> tokens = split_tokens(file.get_field("log_" + tag));
    if (tokens.size() != 3) {
        throw EcoevolityError("checkpoint: malformed log record for " + tag);
    }
    unsigned long long size = string_to_unsigned(tokens.at(0));
    unsigned long long expected_tail_length = string_to_unsigned(tokens.at(1));
    std::string pre = "Cannot resume: the " + tag + " log file \'" + path + "\' ";
    struct stat info;
    if (::stat(path.c_str(), &info) != 0) {
        throw EcoevolityError(pre + "does not exist");
    }
    if (file_size(path) < size) {
        std::ostringstream message;
        message << pre << "is shorter (" << file_size(path)
                << " bytes) than it was when the checkpoint was written ("
                << size << " bytes)";
        throw EcoevolityError(message.str());
    }
    unsigned long tail_length = 0;
    if (hash_file_tail(path, size, tail_length) != tokens.at(2) ||
            (tail_length != expected_tail_length)) {
        throw EcoevolityError(
                pre + "does not match the checkpoint; it was changed after the checkpoint was written");
    }
    return size;
}

// Make sure a log file is the one the checkpoint was taken from, and cut off
// everything after the checkpoint. To leave all the logs alone if any of
// them does not match, check them all with check_log_file first.
inline void roll_back_log_file(
        const CheckpointFile & file,
        const std::string & tag,
        const std::string & path) {
    unsigned long long size = check_log_file(file, tag, path);
    if (file_size(path) > size) {
        truncate_file(path, size);
    }
}

//////////////////////////////////////////////////////////////////////////////
// Everything the MCMC loop needs to know about checkpointing a run
//////////////////////////////////////////////////////////////////////////////

struct Control {
    // Checkpointing is only active if this is true; if it is false, nothing
    // in this struct is used and the MCMC runs exactly as it always has.
    bool enabled = false;
    // Write a checkpoint every this many generations (0: only when stopped by
    // a signal)
    unsigned long every = 0;
    // The checkpoint file
    std::string path;
    // Log files written by the run, which are rolled back to the checkpoint
    // when resuming
    std::string tree_log_path;
    std::string state_log_path;
    std::string operator_log_path;
    // Settings that identify the run, recorded in every checkpoint (and
    // verified when resuming): anything that has to be the same for a
    // continued run to be the run that was stopped
    std::vector< std::pair<std::string, std::string> > identity;
    // Set to continue the run recorded in resume_file
    bool resuming = false;
    CheckpointFile resume_file;
    // Set by the MCMC if it was stopped by a signal
    int stopped_by_signal = 0;
    unsigned long stopped_at_generation = 0;
};


//////////////////////////////////////////////////////////////////////////////
// Graceful stop on SIGINT/SIGTERM (and SIGHUP)
//
// Handlers are only installed when checkpointing is turned on. The first
// signal sets a flag that the MCMC loop polls once per generation. A second
// signal exits immediately. A signal that was being ignored when the run
// started (e.g., SIGINT for a job started in the background by a script, or
// SIGHUP under "nohup") stays ignored.
//
// The flag is a lock-free atomic rather than a volatile sig_atomic_t because,
// with multiple threads, the handler can run on a different thread from the
// one that polls the flag.
//////////////////////////////////////////////////////////////////////////////

// Only lock-free atomics may be used in a signal handler
static_assert(ATOMIC_INT_LOCK_FREE == 2,
        "checkpoint: std::atomic<int> has to be lock-free to be used in a signal handler");

inline std::atomic<int> & signal_flag() {
    static std::atomic<int> flag(0);
    return flag;
}

inline void stop_signal_handler(int signal_number) {
    if (signal_flag().exchange(signal_number) != 0) {
        const char message[] =
                "\nSecond signal received; exiting immediately without writing a new checkpoint.\n";
        ssize_t ignored = ::write(2, message, sizeof(message) - 1);
        (void)ignored;
        ::_exit(EXIT_STATUS_FORCED);
    }
}

struct SavedSignalAction {
    struct sigaction action;
    bool installed;
};

struct SavedSignalActions {
    SavedSignalAction interrupt;
    SavedSignalAction terminate;
    SavedSignalAction hangup;
    bool installed;
};

inline SavedSignalActions & saved_signal_actions() {
    static SavedSignalActions saved;
    return saved;
}

// Install the stop handler for a signal, unless the signal is being ignored
inline void install_stop_signal_handler(
        int signal_number,
        const struct sigaction & action,
        SavedSignalAction & saved) {
    saved.installed = false;
    struct sigaction current;
    if ((sigaction(signal_number, NULL, &current) == 0) &&
            (current.sa_handler != SIG_IGN) &&
            (sigaction(signal_number, &action, &saved.action) == 0)) {
        saved.installed = true;
    }
}

inline void restore_stop_signal_handler(
        int signal_number,
        SavedSignalAction & saved) {
    if (saved.installed) {
        sigaction(signal_number, &saved.action, NULL);
        saved.installed = false;
    }
}

inline void install_stop_signal_handlers() {
    SavedSignalActions & saved = saved_signal_actions();
    if (saved.installed) {
        return;
    }
    signal_flag().store(0);
    struct sigaction action;
    std::memset(&action, 0, sizeof(action));
    action.sa_handler = stop_signal_handler;
    sigemptyset(&action.sa_mask);
    action.sa_flags = 0;
    install_stop_signal_handler(SIGINT, action, saved.interrupt);
    install_stop_signal_handler(SIGTERM, action, saved.terminate);
    install_stop_signal_handler(SIGHUP, action, saved.hangup);
    saved.installed = true;
}

inline void restore_stop_signal_handlers() {
    SavedSignalActions & saved = saved_signal_actions();
    if (! saved.installed) {
        return;
    }
    restore_stop_signal_handler(SIGINT, saved.interrupt);
    restore_stop_signal_handler(SIGTERM, saved.terminate);
    restore_stop_signal_handler(SIGHUP, saved.hangup);
    saved.installed = false;
    signal_flag().store(0);
}

// Zero if no stop was requested; otherwise the signal number
inline int pending_stop_signal() {
    return signal_flag().load();
}

// Intercept stop signals for as long as the object exists (if active)
class StopSignalGuard {
    protected:
        bool active_;

    public:
        explicit StopSignalGuard(bool active) : active_(active) {
            if (this->active_) {
                install_stop_signal_handlers();
            }
        }
        ~StopSignalGuard() {
            if (this->active_) {
                restore_stop_signal_handlers();
            }
        }
        StopSignalGuard(const StopSignalGuard &) = delete;
        StopSignalGuard & operator=(const StopSignalGuard &) = delete;
};

} // namespace checkpoint

#endif
