import cffi

ffibuilder = cffi.FFI()

ffibuilder.cdef(
    """void calculate_gdf(double **, int, double **, int, int, double, double *);"""
)
ffibuilder.set_source(
    "spinner.simple_nn.utils._libgdf",
    '#include "gdf.h"',
    sources=["spinner/simple_nn/utils/gdf.cpp"],
    source_extension=".cpp",
    include_dirs=["spinner/simple_nn/utils/"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
