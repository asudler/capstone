#ifndef MISCTOOLS_H
#define MISCTOOLS_H

#include <string>
#include <vector>

std::vector<std::vector<double>> read
(
	std::string fname,
    char delimiter=',',
	int skip_header=0
);

#endif // MISCTOOLS_H
