#pragma once

#include "cpp11/R.hpp"
#include <stdexcept>
#include <string>

// Keep these development checks active even in builds defining NDEBUG.
#define NATIVE_CHECK(condition) do { \
    if (!(condition)) { \
        throw std::runtime_error(std::string(__FILE__) + ":" + \
            std::to_string(__LINE__) + ": " #condition); \
    } \
} while (false)
