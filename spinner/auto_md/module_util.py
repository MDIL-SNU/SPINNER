import os, datetime

def move_to_directory(directory):
    os.chdir(directory)

def write_log(log_path, text):
    with open(log_path, 'a') as s:
        s.write(text+'\n')

def write_log_with_timestamp(log_path, text):
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
    msg = f"{timestamp} {text}"
    with open(log_path, 'a') as s:
        s.write(msg+'\n')
    return timestamp

def create_and_move_to_directory(directory_name):
    os.makedirs(directory_name, exist_ok=True)
    os.chdir(directory_name)

def calculate_elapsed_time(start_timestamp, end_timestamp):
    start_time = datetime.datetime.strptime(start_timestamp, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end_timestamp, "%Y-%m-%d %H:%M:%S")
    elapsed_time = end_time - start_time
    return elapsed_time

def calculate_total_elapsed_time(log):
    total_time = datetime.timedelta()
    with open(log, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'time:' in line:
                time_to_add = line.split()[-1]
                try:
                    hours, minutes, seconds = map(int, time_to_add.split(':'))
                    elapsed_time = datetime.timedelta(hours=hours, minutes=minutes, seconds=seconds)
                    total_time += elapsed_time
                except:
                    pass
    return total_time

def get_abs_path(*paths):
    path = os.path.join(*paths)
    if not path.startswith('/'):
        path = '/' + path
    return path

def get_divisor(num):
    divisors = []
    for i in range(1, num+1):
        if num%i == 0:
            divisors.append(i)
    return divisors

def get_cpu_core_number():
    return len(os.sched_getaffinity(0))


if __name__ == "__main__":
    print(get_divisor(get_cpu_core_number()))
