import os,sys

from distutils.core import setup, Extension

sys.path.append("/home/kang1717/.local/lib/python3.8/site-packages/pybind11/include")

setup(name='pyrandspg',
      description='randSpg python',
      ext_modules=[Extension("pyrandspg",["pyrandspg.cpp"])]
)
