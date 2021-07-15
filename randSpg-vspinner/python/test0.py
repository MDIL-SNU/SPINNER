import pyrandspg
input_ = pyrandspg.RandSpgInput(3, [4,2,2], pyrandspg.LatticeStruct(1.0, 1.0, 1.0, 60.0, 60.0, 60.0), pyrandspg.LatticeStruct(4.0, 4.0, 4.0, 120.0, 120.0, 120.0), 1.0,  1000.0, 100, 0.0)

c = pyrandspg.RandSpg.randSpgCrystal(input_)

print (c.getPOSCARString())
