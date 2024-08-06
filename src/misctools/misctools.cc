#include <cstring>
#include <iostream>
#include <stdexcept>
#include "misctools.h"

std::vector<std::vector<double>> csv_read::read(std::string fname)
{
    FILE *fp = fopen(fname.c_str(), "r");
    if(fp == NULL) throw std::invalid_argument("file not found");

    char buffer[1024];
    std::vector<std::vector<double>> v;
    while(fgets(buffer, sizeof(buffer), fp))
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
    fclose(fp);

    return v; 
} // csv_read

