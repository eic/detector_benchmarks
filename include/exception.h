#ifndef UTIL_EXCEPTION_H
#define UTIL_EXCEPTION_H

#include <exception>
#include <string>

namespace eic::util {
  class Exception : public std::exception {
  public:
    Exception(std::string_view msg, std::string_view type = "exception") : msg_{msg}, type_{type} {}

    virtual const char* what() const throw() { return msg_.c_str(); }
    virtual const char* type() const throw() { return type_.c_str(); }
    virtual ~Exception() throw() {}

  private:
    std::string msg_;
    std::string type_;
  };
} // namespace eic::util

#endif
