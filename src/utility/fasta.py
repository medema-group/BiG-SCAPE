import os

def fasta_parser(handle):
    """Parses a fasta file, and stores it in a dictionary.
    Only works if there are no duplicate fasta headers in the fasta file"""
    
    fasta_dict = {}
    header = ""
    for line in handle:
        if line[0] == ">":
            header=line.strip()[1:]
        else:
            try:
                fasta_dict[header] += line.strip()
            except KeyError:
                fasta_dict[header] = line.strip()

    return fasta_dict


def get_fasta_keys(handle):
    """Parses a fasta file, only stores headers
    """
    
    header_list = []
    for line in handle:
        if line[0] == ">":
            header_list.append(line.strip()[1:])
            
    return header_list

def save_domain_seqs(filtered_matrix, fasta_dict, domains_folder, outputbase):
    """Write fasta sequences for the domains in the right pfam-domain file"""
    for row in filtered_matrix:
        domain = row[5]
        header = row[-1].strip()
        seq = fasta_dict[header] #access the sequence by using the header
        
        domain_file = open(os.path.join(domains_folder, domain + ".fasta"), 'a') #append to existing file
        domain_file.write(">{}:{}:{}\n{}\n".format(header, row[3], row[4],
            seq[int(row[3]):int(row[4])])) #only use the range of the pfam domain within the sequence
        domain_file.close()
