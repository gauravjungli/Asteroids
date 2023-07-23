#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

std::string doubleToStringWithPrecision(double value, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

int main() {
    double num = 3.141592653589793;
    int precision = 15;

    std::string result = doubleToStringWithPrecision(num, precision);
    std::cout << "Result: " << result << std::endl;

    return 0;
}
