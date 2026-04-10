#pragma once
#include <exception>
#ifndef DISABLE_SPDLOG
#include <spdlog/spdlog.h>
#else
#include <iostream>
#include <sstream>
#endif
#include <string_view>
namespace coacd
{
    namespace logger
    {

        #ifndef DISABLE_SPDLOG
        std::shared_ptr<spdlog::logger> get();
        #else
        // Minimal {}-style format: replaces each {} or {:.Xf} placeholder
        // with the next argument, printed via operator<<.
        inline void fmtApply(std::ostream& out, std::string_view fmt) {
            out << fmt;
        }
        template <typename Arg, typename... Args>
        void fmtApply(std::ostream& out, std::string_view fmt, const Arg& arg, const Args&... rest) {
            auto pos = fmt.find('{');
            if (pos == std::string_view::npos) { out << fmt; return; }
            out << fmt.substr(0, pos);
            // Skip past the closing '}'
            auto end = fmt.find('}', pos);
            if (end == std::string_view::npos) { out << arg; return; }
            out << arg;
            fmtApply(out, fmt.substr(end + 1), rest...);
        }
        template <typename... Args>
        void doPrint(std::ostream& out, std::string_view fmt, const Args&... args) {
            fmtApply(out, fmt, args...);
            out << std::endl;
        }
        #endif

        template <typename... Args>
        inline void debug(std::string_view fmt, const Args &...args)
        {
            #ifndef DISABLE_SPDLOG
            get()->debug(fmt, args...);
            #else
            doPrint(std::cout, fmt, args...);
            #endif
        };

        template <typename... Args>
        inline void info(std::string_view fmt, const Args &...args)
        {
            #ifndef DISABLE_SPDLOG
            get()->info(fmt, args...);
            #else
            doPrint(std::cout, fmt, args...);
            #endif
        };

        template <typename... Args>
        inline void warn(std::string_view fmt, const Args &...args)
        {
            #ifndef DISABLE_SPDLOG
            get()->warn(fmt, args...);
            #else
            doPrint(std::cout, fmt, args...);
            #endif
        };

        template <typename... Args>
        inline void error(std::string_view fmt, const Args &...args)
        {
            #ifndef DISABLE_SPDLOG
            get()->error(fmt, args...);
            #else
            doPrint(std::cout, fmt, args...);
            #endif
        };

        template <typename... Args>
        inline void critical(std::string_view fmt, const Args &...args)
        {
            #ifndef DISABLE_SPDLOG
            get()->critical(fmt, args...);
            #else
            doPrint(std::cout, fmt, args...);
            #endif
        };

    } // namespace logger
} // namespace coacd
