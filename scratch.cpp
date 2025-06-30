#include <iostream>
#include <fstream>
#include <string>

void printMemoryUsage(const std::string& label = "") {
    std::ifstream stat_stream("/proc/self/status", std::ios_base::in);
    std::string line;
    while (std::getline(stat_stream, line)) {
        if (line.find("VmRSS:") == 0 || line.find("VmSize:") == 0) {
            std::cout << label << line << std::endl;
        }
    }
}




