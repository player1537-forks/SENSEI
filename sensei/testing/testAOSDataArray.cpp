#include "AOSDataArray.h"
#include "AOSDataArray.txx"


using sensei::DataArray;
using allocator = DataArray::allocator;

int main(int argc, char **argv)
{
    sensei::AOSDataArray<double> *da = sensei::AOSDataArray<double>::New(allocator::malloc);
    da->Delete();

    return 0;
}
