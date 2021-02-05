/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#pragma once

using clo_arg_t = const char *;
using clo_argv_t = const clo_arg_t *;
using clo_value_t = const char *;

class CmdLineOptions {
public:
    /**
    map a long human readable option name to a single character option.
    the character option must be specified in the options string.
    CmdLineOptions::LongFormat cmd_line_options[] = {
        { '?', "help" },
        { 'a', "alice" },
        { 'b', "bob" },
        { 0, nullptr },
    usage like this: --help --alice required --bob optional
    **/
    class LongFormat {
    public:
        char option_;
        const char *long_name_;
    };

    CmdLineOptions(int argc, clo_argv_t argv, const char *options,
        const LongFormat *lf = nullptr);
    ~CmdLineOptions() noexcept;

    /**
    walk the list of arguments looking for properly formatted arguments.
    this function is similar to gnu's getopt function.
    options are single letters.
    options can be separated like this: -a -b -c -d
    options can be combined like this: -abcd
    -- stops option parsing
    options followed by : require a value
      foo -a hello
      foo -aworld
    options followed by :: may have a value
    parsing stops at an unknown option.
    parsing stops if a required option is missing.
    arg_index_ is the index of the first non-option argument.
    **/
    bool get();

    /** reset the internal state so we can parse the command line again. **/
    void reset();

    char option_;
    clo_value_t value_;
    int arg_index_;
    bool error_;

private:
    int argc_;
    clo_argv_t argv_;
    const char *options_;
    const LongFormat *long_format_;
    clo_arg_t place_;
};
