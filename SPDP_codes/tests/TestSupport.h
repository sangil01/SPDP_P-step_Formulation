#ifndef SPDP_TEST_SUPPORT_H
#define SPDP_TEST_SUPPORT_H

#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace spdp::test {

struct TestCase {
    std::string name;
    std::function<void()> run;
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> cases;
    return cases;
}

class Registrar {
public:
    Registrar(std::string name, std::function<void()> run) {
        registry().push_back(TestCase{std::move(name), std::move(run)});
    }
};

template <typename Left, typename Right>
void check_equal(
    const Left& left,
    const Right& right,
    const char* left_expression,
    const char* right_expression,
    const char* file,
    int line
) {
    if (left == right) {
        return;
    }
    std::ostringstream message;
    message << file << ':' << line << ": expected " << left_expression
            << " == " << right_expression;
    throw std::runtime_error(message.str());
}

inline void check_true(
    bool condition,
    const char* expression,
    const char* file,
    int line
) {
    if (condition) {
        return;
    }
    std::ostringstream message;
    message << file << ':' << line << ": expected " << expression;
    throw std::runtime_error(message.str());
}

}  // namespace spdp::test

#define SPDP_TEST(name)                                                       \
    void name();                                                              \
    static ::spdp::test::Registrar name##_registrar(#name, name);            \
    void name()

#define SPDP_CHECK(expression)                                                \
    ::spdp::test::check_true(                                                 \
        static_cast<bool>(expression), #expression, __FILE__, __LINE__       \
    )

#define SPDP_CHECK_EQ(left, right)                                            \
    ::spdp::test::check_equal(                                                \
        (left), (right), #left, #right, __FILE__, __LINE__                   \
    )

#endif
