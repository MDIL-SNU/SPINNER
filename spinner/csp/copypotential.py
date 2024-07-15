import shutil
import sys

workdirectory = sys.argv[1]
input_dir     = sys.argv[2]
num_of_iteration = int(sys.argv[3])

shutil.copy("potential_saved_iteration"+str(num_of_iteration+1),workdirectory+"/"+input_dir+"/potential")
