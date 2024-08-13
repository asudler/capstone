#include <cstring>
#include <iostream>
#include <stdexcept>
#include "misctools.h"

std::vector<std::vector<double>> read(
	std::string fname,
	int skip_header
)
{
    FILE *fp = fopen(fname.c_str(), "r");
    if(fp == NULL) throw std::invalid_argument("file not found");

    char buffer[1024];
    std::vector<std::vector<double>> v; int counter = 0;
    while(fgets(buffer, sizeof(buffer), fp))
    {
        if(counter > skip_header)
        {
            char *data = std::strtok(buffer, ",");
            std::vector<double> vrow;
            while(data != NULL)
            {
                vrow.push_back(std::stod(data));
                data = strtok(NULL, ",");
            }
        v.push_back(vrow);
        }
        counter += 1;
    } // there is probably a better way to do it... but it works
    fclose(fp);

    return v; 
} // csv_read

