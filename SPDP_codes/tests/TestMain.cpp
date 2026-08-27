#include "TestSupport.h"

#include <exception>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    std::string selected_case;
    if (argc == 3 && std::string(argv[1]) == "--case") {
        selected_case = argv[2];
    } else if (argc != 1) {
        std::cerr << "Usage: SPDP_cp_tests [--case NAME]\n";
        return 2;
    }

    int executed = 0;
    int failed = 0;
    for (const spdp::test::TestCase& test_case : spdp::test::registry()) {
        if (!selected_case.empty() && test_case.name != selected_case) {
            continue;
        }
        ++executed;
        try {
            test_case.run();
            std::cout << "PASS " << test_case.name << '\n';
        } catch (const std::exception& error) {
            ++failed;
            std::cerr << "FAIL " << test_case.name << ": " << error.what()
                      << '\n';
        }
    }

    if (executed == 0) {
        std::cerr << "No matching test case: " << selected_case << '\n';
        return 2;
    }
    return failed == 0 ? 0 : 1;
}
