#ifndef MISCTOOLS_H
#define MISCTOOLS_H

#include <string>
#include <vector>

class csv_read
{
    public:
        csv_read() {}
        virtual ~csv_read() {}
        
        // methods
        static std::vector<std::vector<double>> read(std::string fname);
}; // csv_read

#endif // MISCTOOLS_H
