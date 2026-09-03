#ifndef UTIL_EXCEPTION_H
#define UTIL_EXCEPTION_H

#include <exception>
#include <string>
#include <string_view>

namespace common_bench {
  class Exception : public std::exception {
  public:
    Exception(std::string_view msg, std::string_view type = "exception") : msg_{msg}, type_{type} {}

    const char* what() const noexcept override { return msg_.c_str(); }
    const char* type() const noexcept { return type_.c_str(); }
    ~Exception() noexcept override {}

  private:
    std::string msg_;
    std::string type_;
  };
} // namespace common_bench

#endif
