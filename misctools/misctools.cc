/*#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include "misctools.h"

std::vector<std::vector<double>> read
(
    std::string fname,
    char delimiter,
    int skip_header
)
{
    std::ifstream istrm(fname, std::ios::binary);
    if(!istrm.is_open())
        throw std::invalid_argument("failed to open file");

    std::string line; // skip header lines below
    for(int i = 0; i < skip_header; i++) std::getline(istrm, line);

    std::vector<std::vector<double>> data;
    while(getline(istrm, line))
    {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string word;
        
        while(getline(ss, word, delimiter))
        {
            try
            {
                double number = stod(word); // cast string to double
                row.push_back(number); // add double to row
            }
            catch(const std::invalid_argument& e)
            {
                std::cerr << "failed double cast: " << e.what() << '\n';
            }
        }
        data.push_back(row); // add row to data
    }
    istrm.close();

    return data;
} // read*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include "misctools.h"

std::vector<std::vector<double>> read
(
    std::string fname,
    char delimiter,
    int skip_header
)
{
    std::ifstream istrm(fname, std::ios::binary);
    if(!istrm.is_open())
        throw std::invalid_argument("failed to open file");

    std::string line; // skip header lines below
    for(int i = 0; i < skip_header; i++) std::getline(istrm, line);

    std::vector<std::vector<double>> data;
    while(std::getline(istrm, line))
    {
        // Check if the line is blank (ignoring leading/trailing whitespaces)
        if (line.empty() || 
            line.find_first_not_of(" \t\n\r\f\v") == std::string::npos)
        {
            std::cerr << "Skipping blank line.\n";  
            continue;  
        }

        std::vector<double> row;
        std::stringstream ss(line);
        std::string word;

        while(std::getline(ss, word, delimiter))
        {
            try
            {
                double number = std::stod(word); // cast string to double
                row.push_back(number); // add double to row
            }
            catch(const std::invalid_argument& e)
            {
                std::cerr << "failed double cast: " << e.what() << '\n';
            }
            catch(const std::out_of_range& e)
            {
                std::cerr << "Warning: value " << word << " out of range. "
                          << "Replacing value with 0.\n";
                row.push_back(0.);
            }
        }
        data.push_back(row); // add row to data
    }
    istrm.close();

    return data;
} // read
