import invoke
#from invoke import task

#invoke.run( " g++ $(python-config --cflags) -I/home/kang1717/.local/lib/python3.8/site-packages/pybind11/include -std=c++11 pyrandspg.cpp "
#            "-o libpyrandspg.so "
#)

#invoke.run( " g++  -O3 -Wall -shared -std=c++11 -fPIC "  
#            " -I/home/kang1717/anaconda3/include/python3.8  -I/home/kang1717/.local/lib/python3.8/site-packages/pybind11/include pyrandspg.cpp "
#            "-o libpyrandspg.so"
#)


invoke.run( " g++  -O3 -Wall -shared -std=c++11 -fPIC "  
            "$(python3 -m pybind11 --includes) "
            " -I/home/kang1717/anaconda3/include/python3.8  -I/home/kang1717/.local/lib/python3.8/site-packages/pybind11/include pyrandspg.cpp "
            "-o libpyrandspg.so"
            "`python3-config --extension-suffix` " 
            "-L. -lpyrandspg  -Wl,-rpath,."
)

