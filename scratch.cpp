#include <iostream>
#include <fstream>
#include <string>

void printCurrentMemoryUsageGB() {
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.find("VmRSS:") == 0) {
            std::istringstream iss(line);
            std::string label, value, unit;
            iss >> label >> value >> unit;
            double mem_kb = std::stod(value);
            std::cout << "Current resident memory (RSS / btop MemB): " 
                      << mem_kb / (1024 * 1024) << " GB" << std::endl;
            break;
        }
    }
}




