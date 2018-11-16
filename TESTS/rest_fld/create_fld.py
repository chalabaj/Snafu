import os, re
import shutil
N=0
while os.path.isdir(os.path.join(os.getcwd(),"PREV_RUN"+str(N))):
    print("{} backup folder already exists".format("PREV_RUN"+str(N)))
    N += 1
else:
    backup_folder = os.path.join(os.getcwd(),"PREV_RUN"+str(N))
    os.mkdir(backup_folder)
    print("Creating back-up folder {}".format("PREV_RUN"+str(N)))


rest_files = os.listdir(os.getcwd())
backup_files = []

restart_files = [ rf for rf in rest_files if re.search(r'restart', rf)]
backup_files = backup_files + restart_files
print(backup_files)

for bf in backup_files:
    print(os.path.join(backup_folder,bf))
    shutil.copy(bf,backup_folder)
os.listdir(backup_folder)