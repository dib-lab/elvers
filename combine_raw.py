import os
import os.path
import subprocess
from subprocess import Popen,PIPE

def parse_filename(fastq_filename):

#Astrangia
# LR-Ast4_S5_L004_R2_001.fastq.gz
# 0=sample ID
# 1=sample number
# 2=lane
# 3=read
# 4=file number.fastq.gz

    listoffilestomerge=[]
    fields=fastq_filename.split("_")
    sampleID=fields[0]
    sampleNum=fields[1]
    lane=fields[2]
    read=fields[3]
    extension=fields[4]
    sample=(sampleID,sampleNum,read)
    return sample

def get_newfilename(sample,out_dir):
	newfilename_base="_".join(sample)
	newfilename=out_dir+newfilename_base+".fastq.gz"
	return newfilename

def collect_files(data_dir):
    filename_dictionary={}
    files = os.listdir(data_dir)
    #for root,dirs,files in os.walk(data_dir):
    for filename in files:
        if filename.endswith(".fastq.gz"):
            sample=parse_filename(filename)
            if sample in filename_dictionary.keys():
                filename_dictionary[sample].append(data_dir+filename)
            else:
                filename_dictionary[sample]=[data_dir+filename]
    return filename_dictionary 	

def get_file_string(fileslist):
    sorted_fileslist=sorted(fileslist)
    files_string=" ".join(sorted_fileslist)
    return files_string

def combine_files(filename_dictionary,output_dir):
    for sample in filename_dictionary.keys():
         newfilename=get_newfilename(sample,output_dir)
         print(newfilename)
         files_string=get_file_string(filename_dictionary[sample])
         combine_string="cat "+files_string+" >> "+newfilename
         print(combine_string)
         s=subprocess.Popen(combine_string, shell=True)
         s.wait()

output_dir="/workspace/rosenlab/astrangia/combined/"
data_dir="/groups/rosenlab/loretta/"
filename_dictionary=collect_files(data_dir)
print(filename_dictionary)
for sample in filename_dictionary.keys():
    print(sample)
    print("This is the number of files to combine:",len(filename_dictionary[sample]))
    print(sorted(filename_dictionary[sample]))
combine_files(filename_dictionary,output_dir)
