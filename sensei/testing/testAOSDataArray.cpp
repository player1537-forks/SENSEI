#include "svtkHAMRDataArray.h"
#include "svtkImageData.h"
#include "svtkPointData.h"

int main(int argc, char **argv)
{
    int nx = 2;
    int ny = 2;
    int nz = 1;

    size_t nTups = nx*ny*nz;
    int nComps = 2;

    svtkImageData *id = svtkImageData::New();
    id->SetDimensions(nx, ny, nz);

    double *ptr = (double*)malloc(nTups*nComps*sizeof(double));
    for (size_t i = 0; i < nTups; ++i)
    {
        for (int j = 0; j < nComps; ++j)
        {
            ptr[nComps*i + j] = nComps*i + j;
        }
    }

    svtkHAMRDataArray<double> *da = svtkHAMRDataArray<double>::New(ptr, nTups,
                                                nComps, Allocator::malloc, -1);

    da->SetName("foo");

    //da->Print(std::cerr);

    id->GetPointData()->AddArray(da);

    da->Delete();

    id->Print(std::cerr);

    id->Delete();

    return 0;
}
