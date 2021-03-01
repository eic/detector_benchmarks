#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "exception.h"
#include <fmt/core.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include <string>

// Bookkeeping of test data to store data of one or more tests in a json file to
// facilitate future accounting.
//
// Usage Example 1 (single test):
// ==============================
// 1. define our test
//      eic::util::Test test1{
//        {{"name", "example_test"},
//        {"title", "Example Test"},
//        {"description", "This is an example of a test definition"},
//        {"quantity", "efficiency"},
//        {"target", "1"}}};
// 2. set pass/fail/error status and return value (in this case .99)
//      test1.pass(0.99)
// 3. write our test data to a json file
//      eic::util::write_test(test1, "test1.json");
//
// Usage Example 2 (multiple tests):
// =================================
// 1. define our tests
//      eic::util::Test test1{
//        {{"name", "example_test"},
//        {"title", "Example Test"},
//        {"description", "This is an example of a test definition"},
//        {"quantity", "efficiency"},
//        {"target", "1"}}};
//      eic::util::Test test2{
//        {{"name", "another_test"},
//        {"title", "Another example Test"},
//        {"description", "This is a second example of a test definition"},
//        {"quantity", "resolution"},
//        {"target", "3."}}};
// 2. set pass/fail/error status and return value (in this case .99)
//      test1.fail(10)
// 3. write our test data to a json file
//      eic::util::write_test({test1, test2}, "test.json");

// Namespace for utility scripts, FIXME this should be part of an independent
// library
namespace eic::util {

  struct TestDefinitionError : Exception {
    TestDefinitionError(std::string_view msg) : Exception(msg, "test_definition_error") {}
  };

  // Wrapper for our test data json, with three methods to set the status
  // after test completion (to pass, fail or error). The default value
  // is error.
  // The following fields should be defined in the definitions json
  // for the test to make sense:
  //  - name: unique identifier for this test
  //  - title: Slightly more verbose identifier for this test
  //  - description: Concise description of what is tested
  //  - quantity: What quantity is tested? Unites of value/target
  //  - target: Target value of <quantity> that we want to reach
  //  - value: Actual value of <quantity>
  //  - weight: Weight for this test (this is defaulted to 1.0 if not specified)
  //  - result: pass/fail/error
  struct Test {
    // note: round braces for the json constructor, as else it will pick the wrong
    //       initializer-list constructur (it will put everything in an array)
    Test(const std::map<std::string, std::string>& definition) : json(definition)
    {
      // std::cout << json.dump() << std::endl;
      // initialize with error (as we don't have a value yet)
      error();
      // Check that all required fields are present
      for (const auto& field :
           {"name", "title", "description", "quantity", "target", "value", "result"}) {
        if (json.find(field) == json.end()) {
          throw TestDefinitionError{
              fmt::format("Error in test definition: field '{}' missing", field)};
        }
      }
      // Default "weight" to 1 if not set
      if (json.find("weight") == json.end()) {
        json["weight"] = 1.0;
      }
    }
    // Set this test to pass/fail/error
    void pass(double value) { update_result("pass", value); }
    void fail(double value) { update_result("fail", value); }
    void error(double value = 0) { update_result("error", value); }

    nlohmann::json json;

  private:
    void update_result(std::string_view status, double value)
    {
      json["result"] = status;
      json["value"]  = value;
    }
  };

  void write_test(const std::vector<Test>& data, const std::string& fname)
  {
    nlohmann::json test;
    for (auto& entry : data) {
      test["tests"].push_back(entry.json);
    }
    std::cout << fmt::format("Writing test data to {}\n", fname);
    std::ofstream output_file(fname);
    output_file << std::setw(4) << test << "\n";
  }
  void write_test(const Test& data, const std::string& fname)
  {
    std::vector<Test> vtd{data};
    write_test(vtd, fname);
  }

} // namespace eic::util

#endif
