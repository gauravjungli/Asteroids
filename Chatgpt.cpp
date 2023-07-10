#include <iostream>
#include <fstream>
#include <vector>

// Function to read a 2D array from a file
std::vector<std::vector<double> > read2DArrayFromFile(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<std::vector<double> > array2D;
    std::vector<double> row;

    double num;
    while (file >> num)
    {   
        row.push_back(num);

        // Check if the row is complete
        if (row.size() == 2) // Change 3 to the desired row size
        {  // std::cout<<num<<"\n";
            array2D.push_back(row);
            row.clear();
        }
    }

    return array2D;
}

int main()
{
    std::string filename = "output/files_20.000000_0.800000/"; // Replace with your file name or path
    filename=filename+"grav.txt";
    // Read the first 2D array from file
    std::vector<std::vector<double> > array1 = read2DArrayFromFile(filename);



    // Print the first array
   std::cout << "Array 1:\n";
    for (const auto& row : array1)
    {
        for (double num : row)
        {
            std::cout << num << " ";
        }
        std::cout << "\n";
    }

    return 0;
}

