#!/usr/bin/env nextflow

/*
How to run:
nextflow -C nucl_to_prot_seq.nf.config run nucl_to_prot_seq.nf
*/

params.contigs='input_fasta'

/* input files */
//contig sequences
contig_files = Channel.fromFilePairs("${params.contigs}/*.fasta",size:1)


process seq_to_6frames_prot{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/6_frame_prot", mode:'link'

  input:
  set sample_id, seq from contig_files
  
  output:
  set sample_id, "${sample_id}_prot.fa" into nucl_seq_to_prot
  
  script:
""" 
gotranseq -s ${seq[0]} -o ${sample_id}_prot.fa -f 6 --trim --numcpu=5
"""
}

process frame_prot_blastp{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blastp", mode:'link'

  input:
  set sample_id, frame_prot from nucl_seq_to_prot

  output:
  set sample_id, frame_prot, "${sample_id}_blastp_fmt6.txt" into blastp_res
  
  script:
""" 
blastp -query ${frame_prot} -db TTV_prot_ORFs -out "${sample_id}_blastp_fmt6.txt" -outfmt 6
"""
}

process selection_blastp_res{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blastp", mode:'link'

  input:
  set sample_id, frame_prot, blastp_fmt6 from blastp_res 
  
  output:
  set sample_id, frame_prot, blastp_fmt6, "${sample_id}_blastp_sel.txt" into blastp_sel
  
  script:
""" 
awk '{ if(\$4 > 500.0) { print }}' ${blastp_fmt6} | grep -i "ORF1" | awk '{print "c"\$0}' | awk 'BEGIN { OFS="\t" } { \$1 = substr( \$1, 2, length(\$1)-2 ); print }' | tr '.' ',' | datamash -g 1 max 3 -f | tr ',' '.' > ${sample_id}_blastp_sel.txt
awk '{ if(\$4 > 100.0) { print }}' ${blastp_fmt6} | grep -i "ORF2" | awk '{print "c"\$0}' | awk 'BEGIN { OFS="\t" } { \$1 = substr( \$1, 2, length(\$1)-2 ); print }' | tr '.' ',' | datamash -g 1 max 3 -f | tr ',' '.' >> ${sample_id}_blastp_sel.txt
awk '{ if(\$4 > 50.0) { print }}' ${blastp_fmt6} | grep -i "ORF3" | awk '{print "c"\$0}' | awk 'BEGIN { OFS="\t" } { \$1 = substr( \$1, 2, length(\$1)-2 ); print }' | tr '.' ',' | datamash -g 1 max 3 -f | tr ',' '.' >> ${sample_id}_blastp_sel.txt
awk '{ if(\$4 > 50.0) { print }}' ${blastp_fmt6} | grep -i "ORF4" | awk '{print "c"\$0}' | awk 'BEGIN { OFS="\t" } { \$1 = substr( \$1, 2, length(\$1)-2 ); print }' | tr '.' ',' | datamash -g 1 max 3 -f | tr ',' '.' >> ${sample_id}_blastp_sel.txt
 
"""
}

process prot_frame_selection_blastp{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blastp", mode:'link'

  input:
  set sample_id, frame_prot, blastp_fmt6, blastp_fmt6_sel from blastp_sel 
  
  output:
  set sample_id, frame_prot, "${sample_id}_prot_frame_selected_blastp.txt" into prot_frame_sel_blastp
  
  script:
""" 
#!/usr/bin/python

f1 = open("$blastp_fmt6_sel","r+")
sel_lines = f1.readlines()

f2 = open("$blastp_fmt6","r+")
res_lines = f2.readlines()

f3 = open("${sample_id}_prot_frame_selected_blastp.txt","a+")

for s_line in sel_lines:
    arr = s_line.split()
    m_arr = arr[1:]
    s_str = ' '.join(map(str, m_arr[:-1]))
    for r_line in res_lines:
        r_arr = r_line.split()
        r_str = ' '.join(map(str, r_arr[1:]))
        if arr[0] in r_arr[0] and s_str == r_str:
            f3.write(r_line)

f1.close()
f2.close()
f3.close()

print("$sample_id")
"""
}

process sort_prot_frame_sel_results{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blastp", mode:'link'

  input:
  set sample_id, frame_prot, prot_sel_blastp from prot_frame_sel_blastp 
  
  output:
  set sample_id, frame_prot, "${sample_id}_blastp_sel_sorted.txt" into prot_fram_sel_sorted
  
  script:
""" 
awk '{ if(\$3 > 40.0) { print }}' ${prot_sel_blastp} | sort -k2 -n > "${sample_id}_blastp_sel_sorted.txt"
"""
}

process delete_pos_overlap{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/blastp", mode:'link'

  input:
  set sample_id, frame_prot, res_sorted from prot_fram_sel_sorted
  
  output:
  set sample_id, frame_prot, "${sample_id}_blastp_sel_sorted_overlap_removed.txt" into overlap_rm
  
  script:
"""
#!/usr/bin/python 

def rm_overlap(sel_lines):
    j = 1
    bool_ = 0
    arr_res = []
    i = 0
    while i < len(sel_lines):
        if i == len(sel_lines)-1:
            j = 0
        else:
            j = i + 1
        x=sel_lines[i].split()
        x_line=sel_lines[i]
        y=sel_lines[j].split()
        y_line=sel_lines[j]
        if x[0] == y[0] and j != 0:
            if x[6] == y[6] and x[7] == y[7]:
                bool_=1
            elif x[6] >= y[6] and x[6] <= y[7]:
                bool_=1
                bool_=0
            elif x[6] <= y[6] and x[6] <= y[7] and x[7] >= y[7]:
                bool_=1
                bool_=0
            elif x[6] <= y[6] and x[7] <= y[7] and x[7] >= y[6]:
                bool_=1
                bool_=0
            elif y[6] >= x[6] and y[6] <= x[7]:
                bool_=1
                bool_=0
            elif y[6] <= x[6] and y[6] <= x[7] and y[7] >= x[7]:
                bool_=1
                bool_=0
            elif y[6] <= x[6] and y[7] <= x[7] and y[7] >= x[6]:
                bool_=1
                bool_=0
            else:
                bool_=0
            if bool_ == 1 and float(x[2]) == float(y[2]):
                print("True")
                if float(x[10]) < float(y[10]):
                    arr_res.append(x_line)
                    i = j
                else:
                    arr_res.append(y_line)
                    i = j
            elif bool_ == 1 and float(x[2]) != float(y[2]):
                if float(x[2]) > float(y[2]):
                    arr_res.append(x_line)
                    i = j
                else:
                    arr_res.append(y_line)
                    i = j
            else:
                arr_res.append(x_line)
        elif x[0] != y[0] and j != 0:
            arr_res.append(x_line)
        else:
            arr_res.append(x_line)
        i += 1
    return arr_res

def main():
    f1 = open("$res_sorted","r+")
    sel_lines = f1.readlines()
    f2 = open("${sample_id}_blastp_sel_sorted_overlap_removed.txt","a+")
    tmp = rm_overlap(sel_lines)

    for i in range(1, 20):
        tmp2 = rm_overlap(tmp)
        tmp = tmp2

    for t in tmp:
        f2.write(t)
    
    f1.close()
    f2.close()

main()
"""
}

process extract_frame_prot_seq{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/6_frame_prot", mode:'link'

  input:
  set sample_id, frame_prot, blastp_r from overlap_rm 
  
  output:
  set sample_id, "${sample_id}_extracted_prot.fa", blastp_r  into ext_prot_seq
  
  script:
""" 
cat ${blastp_r} | awk '{print \$1}' > list_prot_tmp.txt

seqtk subseq ${frame_prot} list_prot_tmp.txt > ${sample_id}_extracted_prot.fa
"""
}

process cut_pos_multiORFs_using_extProt{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/multiORF_fasta", mode:'link'

  input:
  set sample_id, ext_prot, blastp_r from ext_prot_seq 
  
  output:
  set sample_id, "${sample_id}_ext_ORFs.fa"  into multiORF
  
  script:
""" 
#!/usr/bin/python 

def split(seq): 
    return [char for char in seq]

def convert(s): 
    new = "" 
    for x in s: 
        new += x    
    return new 

f1 = open("$ext_prot","r+")
prot_seq = f1.readlines()
f2 = open("$blastp_r","r+")
blastp_res = f2.readlines()
f3 = open("${sample_id}_ext_ORFs.fa", "a+")

for i in range(len(blastp_res)):
    for j in range(len(prot_seq)):
        if blastp_res[i].split()[0] in prot_seq[j]:
            tmp=split(prot_seq[j+1])
            start = int(blastp_res[i].split()[6]) - 1
            stop = int(blastp_res[i].split()[7]) - 1
            tmp2 = tmp[start:stop]
            if '\\n' not in tmp2:
                tmp2.append('\\n')
            header = ">" + "${sample_id}" + "_" + blastp_res[i].split()[0]
            if "orf1" in blastp_res[i].split()[1].lower():
                header = header + "_ORF1"
            elif "orf2" in blastp_res[i].split()[1].lower():
                header = header + "_ORF2"
            elif "orf3" in blastp_res[i].split()[1].lower():
                header = header + "_ORF3"
            elif "orf4" in blastp_res[i].split()[1].lower():
                header = header + "_ORF4"
            header = header + " pos:" + blastp_res[i].split()[6] + "-" + blastp_res[i].split()[7] + " alignment_len:" + blastp_res[i].split()[3] + " aa_len:" + str(len(tmp2) - 1) + "\\n"
            f3.write(header)
            f3.write(convert(tmp2))
f1.close()
f2.close()
f3.close()
"""
}

process divide_by_ORFs{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/final_ORFs_results", mode:'link'

  input:
  set sample_id, ext_ORF  from multiORF
  
  output:
  set sample_id, "${sample_id}_ORF1.fa", "${sample_id}_ORF2.fa", "${sample_id}_ORF3.fa", "${sample_id}_ORF4.fa"  into div_ORF
  
  script:
""" 
cat $ext_ORF | grep -i "ORF1" | awk '{print \$1}' | cut -c2- | awk '!x[\$0]++' > tmp1.txt
sleep 1
cat $ext_ORF | grep -i "ORF2" | awk '{print \$1}' | cut -c2- | awk '!x[\$0]++' > tmp2.txt
sleep 1
cat $ext_ORF | grep -i "ORF3" | awk '{print \$1}' | cut -c2- | awk '!x[\$0]++' > tmp3.txt
sleep 1
cat $ext_ORF | grep -i "ORF4" | awk '{print \$1}' | cut -c2- | awk '!x[\$0]++' > tmp4.txt
sleep 1
seqtk subseq ${ext_ORF} tmp1.txt | awk '!x[\$0]++' > ${sample_id}_ORF1.fa
sleep 1
seqtk subseq ${ext_ORF} tmp2.txt | awk '!x[\$0]++' > ${sample_id}_ORF2.fa
sleep 1
seqtk subseq ${ext_ORF} tmp3.txt | awk '!x[\$0]++' > ${sample_id}_ORF3.fa
sleep 1
seqtk subseq ${ext_ORF} tmp4.txt | awk '!x[\$0]++' > ${sample_id}_ORF4.fa
"""
}
