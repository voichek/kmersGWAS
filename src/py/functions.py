from glob import glob
import os
import sys

def get_file_type_in_dir(dir_name, type_suffix):
    fns = glob("%s/*.%s" % (dir_name, type_suffix))
    if len(fns) != 1:
        raise Exception("File count not = 1 %s %d\n%s"  % (type_suffix, len(fns), dir_name))
    return fns[0]

def get_column(fn, index, sep = "\t"):
    return [line.split(sep)[index] for line in file(fn,"r").read().split("\n")[:-1]]

def create_new_relatedness_matrix(fn_in, fn_out, indices):
    def cut_ind(x, ind):
        return [x[ind[i]] for i in range(len(ind))]
    m = [cut_ind(x.split("\t")[:-1],indices) for x in file(fn_in,"r").read().split("\n")[:-1]]
    m = cut_ind(m, indices)
    # convert back to text
    fout = file(fn_out, "w")
    # in the original file all lines ended with "\t" and there was an empty last line
    # I am trying to keep this standard
    fout.write("\n".join(["\t".join(x)+"\t" for x in m]) + "\n")
    fout.close()

def count_running_gemma():
    return int(subprocess.check_output("ps -c | grep gemma | wc -l" , \
            shell=True, stderr = subprocess.PIPE)[:-1])

def get_best_pval(fn):
    return subprocess.check_output("cat %s | tail -n +2 | cut -f9 | awk '$1 < 0.001' | sort -g | head -n1"\
            % fn, shell=True, stderr = subprocess.PIPE)[:-1]

def dir_exist(d):
    if len(d) > 0 and (d[-1] == "/"):
        d = d[:-1]

    return (len(glob(d)) > 0)

def get_file(fn):
    if glob(fn) == 0:
        s = "Couldn't find %s" % fn
        raise Exception(s)
    else:
        return fn

