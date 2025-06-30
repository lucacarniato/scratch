#include <iostream>
#include <fstream>
#include <string>

void printPeakMemoryUsageGB() {
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.find("VmHWM:") == 0) {  // Peak resident set size
            std::istringstream iss(line);
            std::string label, value, unit;
            iss >> label >> value >> unit;
            double mem_kb = std::stod(value);
            std::cout << "Peak memory usage: " << mem_kb / (1024 * 1024) << " GB" << std::endl;
            break;
        }
    }
}





