from glob import glob
import os
import sys
import subprocess
import shlex

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
    m = [cut_ind(x.split("\t"),indices) for x in file(fn_in,"r").read().split("\n")[:-1]]
    m = cut_ind(m, indices)
    # convert back to text
    fout = file(fn_out, "w")
    fout.write("\n".join(["\t".join(x) for x in m]) + "\n")
    fout.close()

def count_running_gemma(g_handles):
    return len([x for x in g_handles if x.poll() == None])

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

def create_dir(dir_path):
    if dir_exist(dir_path):
        err_msg = "directory allready exist %s" % dir_path
        raise Exception(err_msg)
    os.system("mkdir %s" % dir_path)
    return dir_path

def multiple_phenotypes_per_accession(fn):
    return (int(subprocess.check_output("cat %s | tail -n +2 | cut -f1 | sort | uniq -c | awk '$1>1' | wc -l"\
            % fn, shell=True, stderr = subprocess.PIPE)[:-1]) > 0)

def run_and_log_command(cmd, logger):
    logger.writelines("RUN: " + cmd + "\n")
    os.system(cmd)

def run_gemma_cmd(cmd, maximal_to_run, g_handles, logger):
    while count_running_gemma(g_handles) >= maximal_to_run:
        time.sleep(30)
    g_handles.append(subprocess.Popen(shlex.split(cmd)))
#    run_and_log_command(cmd, logger)
    logger.write("RUNG ({}/{}): {}\n".format(count_running_gemma(g_handles), len(g_handles), cmd)) 

def copy_phenotypes(original_file, dest_file,average_script, logger):
    if multiple_phenotypes_per_accession(original_file):
        logger.writelines("Multiple phenotypes found per accessions -> Averaging\n")
        os.system("cat %s | tail -n +2 | awk -f %s > %s" % \
                (original_file, average_script, dest_file))
    else:
        logger.writelines("Unique phenotype per accession, copying phenotype data to directory\n")
        os.system("cp %s %s" % (original_file, dest_file))

def create_full_fam_file(new_fam_file, base_fam, fam_to_expand, n_phenotypes):
    orig_fam_acc = get_column(base_fam, 0, " ")
    cur_fam_acc = get_column(fam_to_expand, 0, "\t")[1:]
    cur_fam_pheno = [get_column(fam_to_expand, i+1, "\t")[1:] for i in range(n_phenotypes)]
    fout = file(new_fam_file, "w")
    for (i, acc_n) in enumerate(orig_fam_acc):
        fout.write("%s %s 0 0 0" % (acc_n, acc_n))
        if acc_n in cur_fam_acc:
            for p_ind in range(n_phenotypes):
                fout.write(" %s" % cur_fam_pheno[p_ind][cur_fam_acc.index(acc_n)])
        else:
            for p_ind in range(n_phenotypes):
                fout.write(" -9")
        fout.write("\n")
    fout.close()

def calc_best_pvals(dir_gemma_output, filename):
    fout = file(filename, "w")
    fns = glob("%s/*.assoc.txt" % dir_gemma_output)
    res = {}
    for fn in fns:
        base_name = fn.split("/")[-1][:-len(".assoc.txt")]
        cur_best_pval = subprocess.check_output(\
                "cat %s | awk 'BEGIN {s=1} {if($9<s) {s = $9}} END {print -log(s)/log(10)}'" % fn \
                ,shell=True, stderr = subprocess.PIPE)[:-1]
        fout.write("%s\t%s\n" % (base_name, cur_best_pval))
        res[base_name] = float(cur_best_pval)
    fout.close()
    return res
        
def get_threshold_from_perm(best_pvals, prefix, n_permutations, p):
    pvals = []
    for i in range(1,n_permutations+1):
        pvals.append(best_pvals[prefix + str(i)])
    pvals.sort(reverse=True)
    return pvals[int(n_permutations * p)-1]
