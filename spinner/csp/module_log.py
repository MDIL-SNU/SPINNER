import time

def write_start_sentence(log_file):
    log_file.write(time.strftime("%c\n",time.localtime(time.time())))
    log_file.write("\n")
    log_file.write("---------------------------------\n")
    log_file.write("     START Genetic Algorithm     \n")
    log_file.write("---------------------------------\n")
    log_file.write("\n")
    log_file.write("\n")

def write_generation_start(log_file,gen):
    log_file.write("----------------   ")
    log_file.write("Generation %d START:  "%(gen))
    log_file.write(time.strftime("%c\n",time.localtime(time.time())))
    log_file.write("\n")
    
def write_generation_end(log_file,gen):
    log_file.write("----------------   ")
    log_file.write("Generation %d END:  "%(gen))
    log_file.write(time.strftime("%c\n",time.localtime(time.time())))
    log_file.write("\n")
