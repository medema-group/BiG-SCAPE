def domtable_parser(gbk, dom_file):
    """Parses the domain table output files from hmmscan"""
    
##example from domain table output:
# target name        accession   tlen query name                                    accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- -----                          -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
#Lycopene_cycl        PF05834.8    378 loc:[0:960](-):gid::pid::loc_tag:['ctg363_1'] -            320   3.1e-38  131.7   0.0   1   1   1.1e-40   1.8e-36  126.0   0.0     7   285    33   295    31   312 0.87 Lycopene cyclase protein
    
#pfd_matrix columns: gbk filename - score - gene id - first coordinate - second coordinate - pfam id - domain name -start coordinate of gene - end coordinate of gene - cds header
    
    pfd_matrix = []
    try:
        dom_handle = open(dom_file, 'r')
    except IOError:
        print("  Error! Could not find file " + dom_file)
    else:
        for line in dom_handle:        
            if line[0] != "#":
                splitline = line.split()
                pfd_row = []
                pfd_row.append(gbk)         #add clustername or gbk filename
                
                pfd_row.append(splitline[13]) #add the score

                header_list = splitline[3].split(":")
                try:
                    pfd_row.append(header_list[header_list.index("gid")+1]) #add gene ID if known
                except ValueError:
                    print("No gene ID in " + gbk)
                    pfd_row.append('')
                    
                pfd_row.append(splitline[19])#first coordinate, env coord from
                pfd_row.append(splitline[20])#second coordinate, env coord to
                
                #===================================================================
                # loc_split = header_list[2].split("]") #second coordinate (of CDS) and the direction
                # pfd_row.append(loc_split[1]) #add direction          
                #===================================================================
                
                pfd_row.append(splitline[1]) #pfam id
                pfd_row.append(splitline[0]) #domain name
                pfd_row.append(header_list[header_list.index("loc")+1]) #start coordinate of gene
                pfd_row.append(header_list[header_list.index("loc")+2]) #end coordinate of gene
                
                pfd_row.append(splitline[3])#cds header
                pfd_matrix.append(pfd_row)

    return pfd_matrix
