/*
Copyright (C) 2012-2021 tim cotter. All rights reserved.
*/

#include "cmd_line.h"

#include <cstring>


CmdLineOptions::CmdLineOptions(
    int argc,
    clo_argv_t argv,
    const char *options,
    const LongFormat *lf
) :
    argc_(argc),
    argv_(argv),
    options_(options),
    long_format_(lf)
{
    reset();
}

CmdLineOptions::~CmdLineOptions() {
}

void CmdLineOptions::reset() {
    option_ = 0;
    value_ = nullptr;
    arg_index_ = 1;
    error_ = false;
    place_ = nullptr;
}

bool CmdLineOptions::get() {
    /** sanity checks. **/
    if (arg_index_ < 0 || argv_ == nullptr || options_ == nullptr) {
        error_ = true;
    }

    /** continue returning errors. **/
    if (error_) {
        option_ = 0;
        value_ = nullptr;
        place_ = nullptr;
        return false;
    }

    /** nothing left to parse. **/
    if (arg_index_ >= argc_) {
        option_ = 0;
        value_ = nullptr;
        place_ = nullptr;
        return false;
    }

    /** attempt to pick up where we left off. **/
    auto cp = place_;
    auto option = '\0';
    if (cp == nullptr || *cp == 0) {
        /**
        no more options in the previous argument.
        look at the next argument.
        **/
        cp = argv_[arg_index_];
        if (cp == nullptr || *cp == 0) {
            /** unexpected end of input. **/
            option_ = 0;
            value_ = nullptr;
            error_ = true;
            place_ = nullptr;
            return false;
        }

        /** option groups start with dash. **/
        option = *cp++;
        if (option != '-') {
            /** end of options. **/
            option_ = 0;
            value_ = nullptr;
            place_ = nullptr;
            return false;
        }

        /** naked dash. **/
        option = *cp++;
        if (option == 0) {
            option_ = 0;
            value_ = nullptr;
            place_ = nullptr;
            return false;
        }

        if (option == '-' && cp[0] == 0) {
            /**
            -- ends parsing.
            arguments begin at the next index.
            **/
            option_ = 0;
            value_ = nullptr;
            ++arg_index_;
            place_ = nullptr;
            return false;
        }

        /** check for long format match. **/
        if (long_format_ && option == '-') {
            for (auto *lf = long_format_; ; ++lf) {
                if (lf->long_name_ == nullptr) {
                    /** long name not in list. **/
                    option_ = 0;
                    value_ = nullptr;
                    error_ = true;
                    place_ = nullptr;
                    return false;
                }
                int len = (int) std::strlen(lf->long_name_);
                if (std::memcmp(cp, lf->long_name_, len) == 0) {
                    cp += len;
                    auto ch = *cp;
                    if (ch == 0) {
                        /** found the long name. **/
                        option = lf->option_;
                        break;
                    }
                    if (ch == '=') {
                        /** found the long name. **/
                        option = lf->option_;
                        break;
                    }
                }
            }
        }
    } else {
        option = *cp++;
    }
    if (option == ':') {
        /** : is not a legal option **/
        option_ = ':';
        value_ = nullptr;
        error_ = true;
        place_ = nullptr;
        return false;
    }

    auto op = strchr(options_, option);
    if (op == nullptr) {
        /** unknown option. **/
        option_ = option;
        value_ = nullptr;
        error_ = true;
        place_ = nullptr;
        return false;
    }

    /** we have found a legal option. **/
    option_ = option;

    /** save our place. **/
    if (*cp) {
        place_ = cp;
    } else {
        ++arg_index_;
        place_ = nullptr;
    }

    /** this option needs no value. **/
    if (op[1] != ':') {
        value_ = nullptr;
        return true;
    }

    /** easy case -c=hello **/
    if (*cp =='=') {
        value_ = cp + 1;
        ++arg_index_;
        place_ = nullptr;
        return true;
    }

    /** easy case -chello **/
    if (*cp) {
        value_ = cp;
        ++arg_index_;
        place_ = nullptr;
        return true;
    }

    /** required value -c hello **/
    if (op[2] != ':') {
        /** no more arguments. **/
        if (arg_index_ >= argc_) {
            value_ = nullptr;
            place_ = nullptr;
            error_ = true;
            return false;
        }

        value_ = argv_[arg_index_];
        if (value_ == nullptr) {
            /** unexpected end of arguments. **/
            error_ = true;
            place_ = nullptr;
            return false;
        }

        ++arg_index_;
        place_ = nullptr;
        return true;
    }

    /** optional value -d **/
    if (arg_index_ >= argc_) {
        value_ = nullptr;
        place_ = nullptr;
        return true;
    }
    cp = argv_[arg_index_];
    if (cp == nullptr) {
        /** unexpected end of arguments. **/
        value_ = nullptr;
        error_ = true;
        place_ = nullptr;
        return true;
    }

    /** no optional value -d -e **/
    if (*cp == '-') {
        value_ = nullptr;
        place_ = nullptr;
        return true;
    }

    /** -d hello **/
    ++arg_index_;
    value_ = cp;
    place_ = nullptr;
    return true;
}
