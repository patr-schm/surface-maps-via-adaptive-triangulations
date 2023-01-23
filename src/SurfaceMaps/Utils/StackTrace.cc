/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Janis Born
 */
#include "StackTrace.hh"

#include <csignal>
#include <iostream>

#if defined(__GLIBCXX__) || defined(__GLIBCPP__)
// GCC: Implement demangling using cxxabi
#include <cxxabi.h>
namespace SurfaceMaps
{
std::string demangle(const std::string& _symbol)
{
    int status;
    char* demangled = abi::__cxa_demangle(_symbol.c_str(), nullptr, nullptr, &status);
    if (demangled)
    {
        std::string result{demangled};
        free(demangled);
        if (status == 0)
            return result;
        else
            return _symbol;
    }
    else
    {
        return _symbol;
    }
}
}
#else
// Other compiler environment: no demangling
namespace SurfaceMaps
{
std::string demangle(const std::string& _symbol)
{
    return _symbol;
}
}
#endif

#ifdef __unix__
#include <execinfo.h>
#include <regex>
namespace SurfaceMaps
{
void print_stack_trace()
{
    void *addresses[20];
    char **strings;

    int size = backtrace(addresses, 20);
    strings = backtrace_symbols(addresses, size);
    std::cerr << "Stack frames: " << size << std::endl;
    // Line format:
    // <path>(<mangled_name>+<offset>) [<address>]
    std::regex line_format{R"(^\s*(.+)\((([^()]+)?\+(0x[0-9a-f]+))?\)\s+\[(0x[0-9a-f]+)\]\s*$)"};
    for(int i = 0; i < size; i++)
    {
        std::string line{strings[i]};
        std::smatch match;
        std::regex_match(line, match, line_format);
        if (!match.empty())
        {
            auto file_name = match[1].str();
            auto symbol = demangle(match[3].str());
            auto offset = match[4].str();
            auto address = match[5].str();
            std::cerr << i << ":";
            if (!file_name.empty()) std::cerr << " " << file_name << " ::";
            if (!symbol.empty()) std::cerr << " " << symbol;
            // Uncomment this if you're interested in address info:
            // if (!offset.empty()) std::cerr << " (+" << offset << ")";
            // if (!address.empty()) std::cerr << " [" << address << "]";
            std::cerr << std::endl;
        }
    }
    free(strings);
}
}
#else
namespace SurfaceMaps
{
void print_stack_trace()
{
    // No-op.
}
}
#endif

namespace SurfaceMaps
{

[[noreturn]]
void handle_segfault(int)
{
    // Prevent infinite recursion by resetting all signals to their default handler (SIG_DFL)
    std::signal(SIGTERM, SIG_DFL);
    std::signal(SIGSEGV, SIG_DFL);
    std::signal(SIGINT,  SIG_DFL);
    std::signal(SIGILL,  SIG_DFL);
    std::signal(SIGABRT, SIG_DFL);
    std::signal(SIGFPE,  SIG_DFL);
    print_stack_trace();
    std::abort();
}

void register_segfault_handler()
{
    std::signal(SIGSEGV, handle_segfault);
    std::signal(SIGILL,  handle_segfault);
    std::signal(SIGABRT, handle_segfault);
    std::signal(SIGFPE,  handle_segfault);
    std::cout << "Registered handler for SIGSEGV, SIGILL, SIGABRT, SIGFPE: handle_segfault" << std::endl;
}

}
