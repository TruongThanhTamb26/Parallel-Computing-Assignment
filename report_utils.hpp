#pragma once

#include <fstream>
#include <iomanip>
#include <string>

inline void appendReport(const std::string& algorithm,
                         const std::string& mode,
                         int rowsA,
                         int colsA,
                         int rowsB,
                         int colsB,
                         double seconds,
                         double checksum,
                         int processes = -1,
                         int threads = -1) {
    std::ofstream report("result_report.txt", std::ios::app);
    if (!report.is_open()) {
        return;
    }

    report << "========================================\n";
    report << "Algorithm : " << algorithm << '\n';
    report << "Mode      : " << mode << '\n';
    report << "Matrix A  : " << rowsA << " x " << colsA << '\n';
    report << "Matrix B  : " << rowsB << " x " << colsB << '\n';

    if (processes > 0) {
        report << "Processes : " << processes << '\n';
    }
    if (threads > 0) {
        report << "Threads   : " << threads << '\n';
    }

    report << std::fixed << std::setprecision(6);
    report << "Time (s)  : " << seconds << '\n';
    report << "Checksum  : " << checksum << '\n';
    report << std::endl;
}
