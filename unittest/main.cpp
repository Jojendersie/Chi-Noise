#include <fstream>
#include <iostream>

void test_generators();
void test_distributions();
void test_sphsampling();
void test_fields();

// Method to store a squared image into .pfm format.
// This is used in the visualization of some test results
bool writePFM(const char* _name, int _size, const float* _data)
{
    std::ofstream file(_name, std::ios::binary);
    if(!file.bad() && !file.fail())
    {
        file.write("Pf\n", sizeof(char) * 3);
        file << _size << " " << _size << "\n";

        file.write("-1.000000\n", sizeof(char) * 10);

        file.write(reinterpret_cast<const char*>(_data), sizeof(float) * _size * _size);

        return true;
    }
    else
    {
        std::cerr << "Error writing hdr image to " << _name;
        return false;
    }
}

int main()
{
    test_distributions();
 //   test_fields();
    test_generators();
  //  test_sphsampling();
    return 0;
}