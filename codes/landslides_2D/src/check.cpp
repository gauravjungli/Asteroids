// Your First C++ Program

#include <iostream>

template <typename T>
int sign(T num) {
    return (num > 0) - (num < 0); // Returns -1 for negative, 1 for positive, 0 for zero
}

int main() {
    std::cout << sign(-10.5) << "\n"; // Output: -1
    std::cout << sign(0.0) << "\n";   // Output: 0
    std::cout << sign(3.14) << "\n";  // Output: 1
}

