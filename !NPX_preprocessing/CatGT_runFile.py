# Written by Chris (Seong Yeol) An
# san11@jh.edu 
# Brown laboratory at Johns Hopkins University
# 5/28/2024
## This file runs CatGT (a common average referencing method developed by Bill Karsh at Janelia)
# CatGT first applies CAR (see car_command) and then outputs also relevant "edge files"

# tools_path: path holding CatGT
# input_path: where the NPX recording file is (you don't need to provide the full path to the file)
# output_path: where CatGT should spit the files to
# probes: number of probes
# run_name: the name of the file
# car_command: we typically use -gblcar but there are also other options - see readme for more info
# filter_command: we skip filtering here because KS 2.5 will apply a filter
# nidq_command: this is an important line that generates event edge files

import os
import subprocess

tools_path = r'D:\NPX_data'
input_path = r'D:\NPX_data\DATA\08312024-NTCA2_g0'
output_path = r'D:\NPX_data\CatGT_OUT'
probes = 0
event_bits = 0
run_name = '08312024-NTCA2_g0'
car_command = '-gblcar'
filter_command = ''
nidq_command = '-xia=0,0,0,0.130,0.110,1500 -xd=0,0,-1,1,0 -xd=0,0,-1,2,0 -xd=0,0,-1,3,0 -xd=0,0,-1,4,0 -xd=0,0,-1,5,0 -xd=0,0,-1,6,0 -xid=0,0,-1,7,0 -xd=0,0,-1,7,0 -xd=0,0,-1,8,0'  

# #run_catgt('06142024-NTBV6_g0',tools_path,input_path,output_path,probes,event_bits,nidq=False,car='global',filter=None,inverted=False)

command = f'runit.bat -dir={input_path} -run={run_name[0:-3]} -g={run_name[-1]},{run_name[-1]} -t=0,0 -ap -prb={probes} {car_command} {filter_command} {nidq_command} -dest={output_path} -pass1_force_ni_ob_bin'
print(f'Running CatGT for {run_name}...')
print(command)
print(os.system(f"cd /d {os.path.join(tools_path, 'CatGT-win')} & {command}"))
print(f'{run_name} done')

# #nidq_command = '-xa=0,0,0,0.2,0.25,0 -xd=0,0,-1,1,0 -xd=0,0,-1,2,0 -xd=0,0,-1,3,0 -xd=0,0,-1,4,0 -xd=0,0,-1,5,0 -xd=0,0,-1,6,0 -xid=0,0,-1,6,0 -xd=0,0,-1,7,0 -xd=0,0,-1,8,0'  


# # ################

command = f'runit.bat -dir={input_path} -run={run_name[0:-3]} -g={run_name[-1]},{run_name[-1]} -t=0,0 -ni -prb={probes} {car_command} {filter_command} {nidq_command} -dest={output_path} -pass1_force_ni_ob_bin'
print(f'Running CatGT for {run_name}...')
print(command)
print(os.system(f"cd /d {os.path.join(tools_path, 'CatGT-win')} & {command}"))
print(f'{run_name} done')




# ### supercat to concatenate
# input_path1 = r'D:\NPX_data\DATA\supercat'
# input_path2 = r'D:\NPX_data\DATA\supercat'
# supercat_output = r'D:\NPX_data\CatGT_OUT'
# ga1 = '08262024-NTCA7_g0'
# ga2 = '08262024-NTCA7_g1'
# command = f'runit.bat -supercat={input_path1,ga1}{input_path2,ga2} -ni -prb=0 {nidq_command} -supercat_trim_edges -dest={supercat_output}'
# print(f'Running CatGT for {run_name}...')
# print(command)
# print(os.system(f"cd /d {os.path.join(tools_path, 'CatGT-win')} & {command}"))
# print(f'{run_name} done')

