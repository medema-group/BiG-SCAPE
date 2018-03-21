#! /usr/bin/python

######################################################################
#                                                                    #
#           PLOT ARROWS FOR GENE CLUSTER GIVEN A GenBank FILE        #
#                           Peter Cimermancic                        #
#                               April 2010                           #
#                heavily modified by Jorge Navarro 2016              #
######################################################################

# Makes sure the script can be used with Python 2 as well as Python 3.
from __future__ import print_function, division
from sys import version_info
if version_info[0]==2:
    range = xrange

import os
import sys
import argparse
from Bio import SeqIO
from random import uniform
from colorsys import hsv_to_rgb
from colorsys import rgb_to_hsv
from math import sin, atan2, pi
from collections import defaultdict

global internal_domain_margin
global gene_contour_thickness
global stripe_thickness
global gene_categories_color

internal_domain_margin = 3
domain_contour_thickness = 1
gene_contour_thickness = 2 # thickness grows outwards
stripe_thickness = 3

domains_color_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "domains_color_file.tsv")


# read various color data
def read_color_genes_file():
    # Try to read already-generated colors for genes
    color_genes = {}
    
    if os.path.isfile(gene_color_file):
        print("  Found file with gene colors")
        with open(gene_color_file, "r") as color_genes_handle:
            for line in color_genes_handle:
                # handle comments and empty lines
                if line[0] != "#" and line.strip():
                    row = line.strip().split("\t")
                    name = row[0]
                    rgb = row[1].split(",")
                    color_genes[name] = [int(rgb[x]) for x in range(3)]
    else:
        print("  Gene color file was not found. A new file will be created")
        with open(gene_color_file, "w") as color_genes_handle:
            color_genes_handle.write("NoName\t255,255,255\n")
        color_genes = {"NoName":[255, 255, 255]}
    
    return color_genes


def read_color_domains_file():
    # Try to read colors for domains
    color_domains = {}
    
    if os.path.isfile(domains_color_file):
        print("  Found file with domains colors")
        with open(domains_color_file, "r") as color_domains_handle:
            for line in color_domains_handle:
                # handle comments and empty lines
                if line[0] != "#" and line.strip():
                    row = line.strip().split("\t")
                    name = row[0]
                    rgb = row[1].split(",")
                    color_domains[name] = [int(rgb[x]) for x in range(3)]
    else:
        print("  Domains colors file was not found. An empty file will be created")
        color_domains_handle = open(domains_color_file, "a+")
        
    return color_domains


# Try to read categories:
def read_pfam_domain_categories():
    pfam_category = {}
    
    if os.path.isfile(pfam_domain_categories):
        print("  Found file with Pfam domain categories")
        with open(pfam_domain_categories, "r") as cat_handle:            
            for line in cat_handle:
                # handle comments and empty lines
                if line[0] != "#" and line.strip():
                    row = line.strip().split("\t")
                    domain = row[1]
                    category = row[0]
                    pfam_category[domain] = category
    else:
        print("  File pfam_domain_categories was NOT found")
                    
    return pfam_category
   

# --- Draw arrow for gene
def draw_arrow(additional_tabs, X, Y, L, l, H, h, strand, color, color_contour, category, gid, domain_list):
    """
    SVG code for arrow:
        - (X,Y) ... upper left (+) or right (-) corner of the arrow
        - L ... arrow length
        - H ... arrow height
        - strand
        - h ... arrow head edge width
        - l ... arrow head length
        - color
        - strand
    the edges are ABCDEFG starting from (X,Y)     
    """

    if strand == '+':
        head_end = L
        if L < l:
            # squeeze arrow if length shorter than head length
            A = [X,Y-h]
            B = [X+L,Y+H/2]
            C = [X,Y+H+h]
            head_start = 0
            points = [A, B, C]
        else:
            A = [X,Y]
            B = [X+L-l,Y]
            C = [X+L-l,Y-h]
            D = [X+L,Y+H/2]
            E = [X+L-l,Y+H+h]
            F = [X+L-l,Y+H]
            G = [X,Y+H]
            head_start = L - l # relative to the start of the gene, not absolute coords.
            points = [A, B, C, D, E, F, G]

    elif strand == '-':
        head_start = 0
        if L < l:
            # squeeze arrow if length shorter than head length
            A = [X,Y+H/2]
            B = [X+L,Y-h]
            C = [X+L,Y+H+h]
            head_end = L
            points = [A, B, C]
        else:
            A = [X+L,Y]
            B = [X+l,Y]
            C = [X+l,Y-h]
            D = [X,Y+H/2]
            E = [X+l,Y+H+h]
            F = [X+l,Y+H]
            G = [X+L,Y+H]
            head_end = l
            points = [A, B, C, D, E, F, G]

    else:
        return ""
    
    head_length = head_end - head_start
    if head_length == 0:
        return ""
    
    points_coords = []
    for point in points:
        points_coords.append(str(int(point[0])) + "," + str(int(point[1])))
    
    arrow = additional_tabs + "\t<g>\n"
    
    # unidentified genes don't have a title and have a darker contour
    if gid != "NoName":
        arrow += additional_tabs + "\t\t<title>" + gid + "</title>\n"
    else:
        color_contour = [50, 50, 50]
        
    arrow += "{}\t\t<polygon class=\"{}\" ".format(additional_tabs, gid)
    arrow += "points=\"{}\" fill=\"rgb({})\" ".format(" ".join(points_coords), ",".join([str(val) for val in color]))
    arrow += "fill-opacity=\"1.0\" stroke=\"rgb({})\" ".format(",".join([str(val) for val in color_contour]))
    arrow += "stroke-width=\"{}\" {} />\n".format(str(gene_contour_thickness), category)
    
    # paint domains. Domains on the tip of the arrow should not have corners sticking
    #  out of them
    for domain in domain_list:
        #[X, L, H, domain_accession, (domain_name, domain_description), color, color_contour]
        dX = domain[0]
        dL = domain[1]
        dH = domain[2]
        dacc = domain[3]
        dname = domain[4][0]
        ddesc = domain[4][1]
        dcolor = domain[5]
        dccolour = domain[6]
        
        arrow += additional_tabs + "\t\t<g>\n"
        arrow += "{}\t\t\t<title>{} ({})\n\"{}\"</title>\n".format(additional_tabs, dname, dacc, ddesc)
        
        if strand == "+":
            # calculate how far from head_start we (the horizontal guide at y=Y+internal_domain_margin)
            #  would crash with the slope
            # Using similar triangles:
            collision_x = head_length * (h + internal_domain_margin)
            collision_x /= (h + H/2.0)
            collision_x = round(collision_x)
            
            # either option for x_margin_offset work
            #m = -float(h + H/2)/(head_length) #slope of right line
            #x_margin_offset = (internal_domain_margin*sqrt(1+m*m))/m
            #x_margin_offset = -(x_margin_offset)
            x_margin_offset = internal_domain_margin/sin(pi - atan2(h+H/2.0,-head_length))

            if (dX + dL) < head_start + collision_x - x_margin_offset:
                arrow += "{}\t\t\t<rect class=\"{}\" x=\"{}\" ".format(additional_tabs, dacc, str(X+dX))
                arrow += "y=\"{}\" stroke-linejoin=\"round\" ".format(str(Y + internal_domain_margin))
                arrow += "width=\"{}\" height=\"{}\" ".format(str(dL), str(dH))
                arrow += "fill=\"rgb({})\" stroke=\"rgb({})\" ".format(",".join([str(val) for val in dcolor]), ",".join([str(val) for val in dccolour]))
                arrow += "stroke-width=\"{}\" opacity=\"0.75\" />\n".format(str(domain_contour_thickness))
            else:
                del points[:]
                
                if dX < head_start + collision_x - x_margin_offset:
                    # add points A and B
                    points.append([X + dX, Y + internal_domain_margin])
                    points.append([X + head_start + collision_x - x_margin_offset, Y + internal_domain_margin])
                    
                else:
                    # add point A'
                    start_y_offset = (h + H/2)*(L - x_margin_offset - dX)
                    start_y_offset /= head_length
                    start_y_offset = int(start_y_offset)
                    points.append([X + dX, int(Y + H/2 - start_y_offset)])
                    
                    
                # handle the rightmost part of the domain
                if dX + dL >= head_end - x_margin_offset: # could happen more easily with the scaling
                    points.append([X + head_end - x_margin_offset, int(Y + H/2)]) # right part is a triangle
                else:
                    # add points C and D
                    end_y_offset = (2*h + H)*(L - x_margin_offset - dX - dL)
                    end_y_offset /= 2*head_length
                    end_y_offset = int(end_y_offset)

                    points.append([X + dX + dL, int(Y + H/2 - end_y_offset)])
                    points.append([X + dX + dL, int(Y + H/2 + end_y_offset)])
            
                # handle lower part
                if dX < head_start + collision_x - x_margin_offset:
                    # add points E and F
                    points.append([X + head_start + collision_x - x_margin_offset, Y + H - internal_domain_margin])
                    points.append([X + dX, Y + H - internal_domain_margin])                    
                else:
                    # add point F'
                    points.append([X + dX, int(Y + H/2 + start_y_offset)])
            
                       
                del points_coords[:]
                for point in points:
                    points_coords.append(str(int(point[0])) + "," + str(int(point[1])))
                    
                arrow += "{}\t\t\t<polygon class=\"{}\" ".format(additional_tabs, dacc)
                arrow += "points=\"{}\" stroke-linejoin=\"round\" ".format(" ".join(points_coords))
                arrow += "width=\"{}\" height=\"{}\" ".format(str(dL), str(dH))
                arrow += "fill=\"rgb({})\" ".format(",".join([str(val) for val in dcolor]))
                arrow += "stroke=\"rgb({})\" ".format(",".join([str(val) for val in dccolour]))
                arrow += "stroke-width=\"{}\" opacity=\"0.75\" />\n".format(str(domain_contour_thickness))
            
        # now check other direction
        else:
            # calculate how far from head_start we (the horizontal guide at y=Y+internal_domain_margin)
            #  would crash with the slope
            # Using similar triangles:
            collision_x = head_length * ((H/2) - internal_domain_margin)
            collision_x /= (h + H/2.0)
            collision_x = round(collision_x)
            
            x_margin_offset = round(internal_domain_margin/sin(atan2(h+H/2.0,head_length)))
            
            # nice, blocky domains
            if dX > collision_x + x_margin_offset:
                arrow += "{}\t\t\t<rect class=\"{}\" ".format(additional_tabs, dacc)
                arrow += "x=\"{}\" y=\"{}\" ".format(str(X+dX), str(Y + internal_domain_margin))
                arrow += "stroke-linejoin=\"round\" width=\"{}\" height=\"{}\" ".format(str(dL), str(dH))
                arrow += "fill=\"rgb({})\" ".format(",".join([str(val) for val in dcolor]))
                arrow += "stroke=\"rgb({})\" ".format(",".join([str(val) for val in dccolour]))
                arrow += "stroke-width=\"{}\" opacity=\"0.75\" />\n".format(str(domain_contour_thickness))
            else:
                del points[:]
                
                # handle lefthand side. Regular point or pointy start?
                if dX >= x_margin_offset:
                    start_y_offset = round((h + H/2)*(dX - x_margin_offset)/head_length)
                    points.append([X + dX, Y + H/2 - start_y_offset])
                else:
                    points.append([X + x_margin_offset, Y + H/2])
                    
                    
                # handle middle/end
                if dX + dL < collision_x + x_margin_offset:
                    if head_length != 0:
                        end_y_offset = round((h + H/2)*(dX + dL - x_margin_offset)/head_length)
                    else:
                        end_y_offset = 0
                    points.append([X + dX + dL, Y + H/2 - end_y_offset])
                    points.append([X + dX + dL, Y + H/2 + end_y_offset])
                else:
                    points.append([X + collision_x + x_margin_offset, Y + internal_domain_margin])
                    points.append([X + dX + dL, Y + internal_domain_margin])
                    points.append([X + dX + dL, Y + internal_domain_margin + dH])
                    points.append([X + collision_x + x_margin_offset, Y + internal_domain_margin + dH])
                    
                # last point, if it's not a pointy domain
                if dX >= x_margin_offset:
                    points.append([X + dX, Y + H/2 + start_y_offset])
                       
                del points_coords[:]
                for point in points:
                    points_coords.append(str(int(point[0])) + "," + str(int(point[1])))
                    
                arrow += "{}\t\t\t<polygon class=\"{}\" ".format(additional_tabs, dacc)
                arrow += "points=\"{}\" stroke-linejoin=\"round\" ".format(" ".join(points_coords))
                arrow += "width=\"{}\" height=\"{}\" ".format(str(dL), str(dH))
                arrow += "fill=\"rgb({})\" ".format(",".join([str(val) for val in dcolor]))
                arrow += "stroke=\"rgb({})\" ".format(",".join([str(val) for val in dccolour]))
                arrow += "stroke-width=\"{}\" opacity=\"0.75\" />\n".format(str(domain_contour_thickness))
        
        arrow += additional_tabs + "\t\t</g>\n"
    
    arrow += additional_tabs + "\t</g>\n"

    return arrow


def draw_line(X,Y,L):
    """
    Draw a line below genes
    """
    
    line = "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" style=\"stroke:rgb(70,70,70); stroke-width:{} \"/>\n".format(str(X), str(Y), str(X+L), str(Y), str(stripe_thickness))
    
    return line


def new_color(gene_or_domain):
    # see https://en.wikipedia.org/wiki/HSL_and_HSV
    # and http://stackoverflow.com/a/1586291
    
    h = uniform(0, 1) # all possible colors

    if gene_or_domain == "gene":
        s = uniform(0.15, 0.3)
        v = uniform(0.98, 1.0)
    elif gene_or_domain == "domain":
        s = uniform(0.5, 0.75) # lower: less saturated
        v = uniform(0.7, 0.9) # lower = darker
    else:
        sys.exit("unknown kind of color. Should be 'gene' or 'domain'")
        
    r, g, b = tuple(int(c * 255) for c in hsv_to_rgb(h, s, v))
    
    return [r, g, b]


def SVG(write_html, outputfile, GenBankFile, BGCname, pfdFile, use_pfd, color_genes, color_domains, pfam_domain_categories, pfam_info, loci, max_width, H=30, h=15, l=30, mX=10, mY=10, scaling=30, absolute_start=0, absolute_end=-1):
    '''
    Create the main SVG document:
        - read pfd file with domain information (pfdFile contains complete path)
        - read GenBank document (GenBankFile contains handle of opened file)
        - record genes, start and stop positions, and strands, and associate domains
        - write the SVG files
    '''
    
    # for colors not found in colors_genes and color_domains, we need to generate them from scratch
    new_color_genes = {}
    new_color_domains = {}
    
    SVG_TEXT = "" # here we keep all the text that will be written down as a file
    
    # check whether we have a corresponding pfd file wih domain annotations
    if use_pfd:
        if not os.path.isfile(pfdFile):
            sys.exit("Error (Arrower): " + pfdFile + " not found")
   

    # --- create SVG header. We have to get max_width first
    # This means that we have to read the gbk file once to know num loci, max_width
    if loci == -1:
        try:
            records = list(SeqIO.parse(GenBankFile), "genbank")
        except:
            sys.exit(" Arrower: error while opening GenBank")
        else:
            loci = len(records)
            max_width = 0
            for record in records:
                if len(record) > max_width:
                    max_width = len(record)
    
        
    if absolute_end < 0: # absolute_end == -1 means "the whole region"
        absolute_end = max_width
    else:
        if (absolute_end - absolute_start) < max_width: # user specified something shorter than full region
            max_width = float(absolute_end - absolute_start)
        else: # user specified something bigger than full region. Cropping to max_width
            absolute_end = max_width
            
    max_width /= scaling
            
    if write_html:
        header = "\t\t<div title=\"" + BGCname + "\">\n"
        additional_tabs = "\t\t\t"
        
        header += "{}<svg width=\"{}\" height=\"{}\">\n".format(additional_tabs, str(max_width + 2*(mX)), str(loci*(2*h + H + 2*mY)))

        addY = loci*(2*h + H + 2*mY)
    else:
        header = "<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" width=\"" + str(max_width + 2*(mX)) + "\" height=\"" + str(loci*(2*h + H + 2*mY)) + "\">\n"
        addY = 0
        
        additional_tabs = "\t"
              
    SVG_TEXT = header
              
    # For info on the color matrix definition: 
    #  https://www.w3.org/TR/SVG11/filters.html#feColorMatrixElement
    # Core Bio: "#DC143C", (220, 20, 60) Dark red
    # Other Bio: 
    #  original: "#DF809D", (223, 128, 157) Pink .87, 0.5, 0.61
    #  alternative: #f4a236, (244,162,54) 0.95, 0.63, 0.21
    # Transporter: "#3F9FBA" (63, 159, 186) Blue
    #  32839a, (50, 131, 154), 0.19, 0.51, 0.6
    # Regulator: "#63BB6D" (99, 187, 109) Green
    #  #127E1B, (18,126,27) 0.07, 0.49, 0.1
    if len(pfam_domain_categories) > 0:
        filters = additional_tabs + "<filter id=\"shadow_CoreBio\" color-interpolation-filters=\"sRGB\" x=\"-65%\" y=\"-25%\" width=\"230%\" height=\"150%\">\n"
        filters += additional_tabs + "\t<feColorMatrix in=\"SourceGraphic\" result=\"matrixOut\" type=\"matrix\" values=\"0 0 0 0 0.85 0 0 0 0 0.08 0 0 0 0 0.23 0 0 0 1 0\" />\n"
        filters += additional_tabs + "\t<feGaussianBlur in=\"matrixOut\" result=\"blurOut\" stdDeviation=\"7\" />\n"
        filters += additional_tabs + "\t<feBlend in=\"SourceGraphic\" in2=\"blurOut\" mode=\"normal\" />\n"
        filters += additional_tabs + "</filter>\n"
        
        filters += additional_tabs + "<filter id=\"shadow_OtherBio\" color-interpolation-filters=\"sRGB\" x=\"-65%\" y=\"-25%\" width=\"230%\" height=\"150%\">\n"
        filters += additional_tabs + "\t<feColorMatrix in=\"SourceGraphic\" result=\"matrixOut\" type=\"matrix\" values=\"0 0 0 0 0.95 0 0 0 0 0.63 0 0 0 0 0.21 0 0 0 1 0\" />\n"
        filters += additional_tabs + "\t<feGaussianBlur in=\"matrixOut\" result=\"blurOut\" stdDeviation=\"7\" />\n"
        filters += additional_tabs + "\t<feBlend in=\"SourceGraphic\" in2=\"blurOut\" mode=\"normal\" />\n"
        filters += additional_tabs + "</filter>\n"
        
        filters += additional_tabs + "<filter id=\"shadow_Transporter\" color-interpolation-filters=\"sRGB\" x=\"-65%\" y=\"-25%\" width=\"230%\" height=\"150%\">\n"
        filters += additional_tabs + "\t<feColorMatrix in=\"SourceGraphic\" result=\"matrixOut\" type=\"matrix\" values=\"0 0 0 0 0.19 0 0 0 0 0.51 0 0 0 0 0.6 0 0 0 1 0\" />\n"
        filters += additional_tabs + "\t<feGaussianBlur in=\"matrixOut\" result=\"blurOut\" stdDeviation=\"7\" />\n"
        filters += additional_tabs + "\t<feBlend in=\"SourceGraphic\" in2=\"blurOut\" mode=\"normal\" />\n"
        filters += additional_tabs + "</filter>\n"
        
        filters += additional_tabs + "<filter id=\"shadow_Regulator\" color-interpolation-filters=\"sRGB\" x=\"-65%\" y=\"-25%\" width=\"230%\" height=\"150%\">\n"
        filters += additional_tabs + "\t<feColorMatrix in=\"SourceGraphic\" result=\"matrixOut\" type=\"matrix\" values=\"0 0 0 0 0.07 0 0 0 0 0.49 0 0 0 0 0.1 0 0 0 1 0\" />\n"
        filters += additional_tabs + "\t<feGaussianBlur in=\"matrixOut\" result=\"blurOut\" stdDeviation=\"7\" />\n"
        filters += additional_tabs + "\t<feBlend in=\"SourceGraphic\" in2=\"blurOut\" mode=\"normal\" />\n"
        filters += additional_tabs + "</filter>\n"
        
        SVG_TEXT += filters

    # --- read in GenBank file

    # handle domains
    if use_pfd:
        identifiers = defaultdict(list)
        with open(pfdFile, "r") as pfd_handle:
            for line in pfd_handle:
                row = line.strip().split("\t")
                
                # use to access to parent's properties
                identifier = row[9].replace("<","").replace(">","")
                # if it's the new version of pfd file, we can take the last part 
                #  to make it equal to the identifiers used in gene_list. Strand
                #  is recorded in parent gene anyway
                if ":strand:+" in identifier:
                    identifier = identifier.replace(":strand:+", "")
                    strand = "+"
                if ":strand:-" in identifier:
                    identifier = identifier.replace(":strand:-", "")
                    strand = "-"
                

                width = 3*(int(row[4]) - int(row[3]))
                            
                if strand == "+":
                    # multiply by 3 because the env. coordinate is in aminoacids, not in bp
                    # This start is relative to the start of the gene
                    start = 3*int(row[3])
                else:
                    loci_start = int(row[7].replace("<","").replace(">",""))
                    loci_end = int(row[8].replace("<","").replace(">",""))
                                    
                    start = loci_end - loci_start - 3*int(row[3]) - width
                
                # geometry
                start = int(start/scaling)
                width = int(width/scaling)

                # accession
                domain_acc = row[5].split(".")[0]
                
                # colors
                try:
                    color = color_domains[domain_acc]
                except KeyError:
                    color = new_color("domain")
                    new_color_domains[domain_acc] = color
                    color_domains[domain_acc] = color
                    pass
                # contour color is a bit darker. We go to h,s,v space for that
                h_, s, v = rgb_to_hsv(float(color[0])/255.0, float(color[1])/255.0, float(color[2])/255.0)
                color_contour = tuple(int(c * 255) for c in hsv_to_rgb(h_, s, 0.8*v))


                # [X, L, H, domain_acc, color, color_contour]
                identifiers[identifier].append([start, width, int(H - 2*internal_domain_margin), domain_acc, pfam_info[domain_acc], color, color_contour])
    
    loci = 0
    feature_counter = 1
    records = list(SeqIO.parse(GenBankFile, "genbank"))
    for seq_record in records:
        add_origin_Y = loci * (2*(h+mY) + H)

        # draw a line that coresponds to cluster size
        ClusterSize = len(seq_record.seq)
        if (absolute_end - absolute_start) < ClusterSize:
            ClusterSize = (absolute_end - absolute_start)
        
        line = draw_line(mX, add_origin_Y + mY + h + H/2, ClusterSize/scaling)
        
        SVG_TEXT += additional_tabs + "<g>\n"
        
        SVG_TEXT += additional_tabs + "\t" + line
        
        # Calculate features for all arrows
        
        for feature in [feature for feature in seq_record.features if feature.location.start >= absolute_start and feature.location.end <= absolute_end]:
            if feature.type == 'CDS':
                # Get name
                try:
                    GeneName = feature.qualifiers['gene'][0]
                    cds_tag = GeneName
                except KeyError:
                    GeneName = 'NoName'
                    cds_tag = ""
                
                if "locus_tag" in feature.qualifiers:
                    cds_tag += " (" + feature.qualifiers["locus_tag"][0] + ")"
                if "product" in feature.qualifiers:
                    cds_tag += "\n" + feature.qualifiers["product"][0]
                    
                # Get color
                color = (255,255,255)
                #try:
                    #color = color_genes[GeneName]
                #except KeyError:
                    #color = new_color("gene")
                    #new_color_genes[GeneName] = color
                    #color_genes[GeneName] = color
                    #pass
                
                color_contour = (0,0,0)
                # change to hsv color palette to lower shade for contour color
                #h_, s, v = rgb_to_hsv(float(color[0])/255.0, float(color[1])/255.0, float(color[2])/255.0)
                #color_contour = tuple(int(c * 255) for c in hsv_to_rgb(h_, s, 0.4*v))
                
                # Get strand
                strand = feature.strand
                if strand == -1:
                    strand = '-'
                elif strand == 1:
                    strand = '+'
                else:
                    sys.exit("Weird strand value: " + strand)
                
                # define arrow's start and end
                # http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html#start
                start = feature.location.start - absolute_start
                start = int(start/scaling)
                stop = feature.location.end - absolute_start
                stop = int(stop/scaling)
                
                # assemble identifier to match domains with this feature
                try:
                    protein_id = feature.qualifiers['protein_id'][0]
                except KeyError:
                    protein_id = ""
                    pass
                identifier = BGCname + "_ORF" + str(feature_counter)
                identifier += ":gid::" if GeneName == "NoName" else ":gid:" + str(GeneName) + ":"
                identifier += "pid:" + str(protein_id) + ":loc:" + str(feature.location.start) + ":" + str(feature.location.end)
                identifier = identifier.replace("<","").replace(">","")

                # gene category according to domain content
                #has_core = False
                #has_otherbio = False
                #has_transporter = False
                #has_regulator = False
                #for row in identifiers[identifier]:
                    #dom_acc = row[3]
                    #cat = ""
                    #try:
                        #cat = pfam_domain_categories[dom_acc]
                    #except KeyError:
                        #pass
                    
                    #if cat == "Core Biosynthetic":
                        #has_core = True
                    #if cat == "Other Biosynthetic":
                        #has_otherbio = True
                    #if cat == "Transporter":
                        #has_transporter = True
                    #if cat == "Regulator":
                        #has_regulator = True
                            
                gene_category = ""
                #if has_core:
                    #gene_category = "filter=\"url(#shadow_CoreBio)\""
                #if has_otherbio and not (has_core or has_transporter or has_regulator):
                    #gene_category = "filter=\"url(#shadow_OtherBio)\""
                #if has_transporter and not (has_core or has_otherbio or has_regulator):
                    #gene_category = "filter=\"url(#shadow_Transporter)\""
                #if has_regulator and not (has_core or has_otherbio or has_transporter):
                    #gene_category = "filter=\"url(#shadow_Regulator)\""
                        
                
                #X, Y, L, l, H, h, strand, color, color_contour, category, gid, domain_list
                arrow = draw_arrow(additional_tabs, start+mX, add_origin_Y+mY+h, int(feature.location.end-feature.location.start)/scaling, l, H, h, strand, color, color_contour, gene_category, cds_tag, identifiers[identifier])
                if arrow == "":
                    print("  (ArrowerSVG) Warning: something went wrong with {}".format(BGCname))
                SVG_TEXT += arrow
                
                feature_counter += 1
                
        loci += 1
        
        SVG_TEXT += additional_tabs + "</g>\n"

    SVG_TEXT += additional_tabs[:-2] + "</svg>\n"
    
    if write_html:
        SVG_TEXT += "\t\t</div>\n"
    
    # finally append new colors to file:
    #if len(new_color_genes) > 0:
        #if len(new_color_genes) < 10:
            #print("  Saving new color names for genes " + ", ".join(new_color_genes.keys()))
        #else:
            #print("  Saving new color names for 10+ genes...")
            
        #with open(gene_color_file, "a") as color_genes_handle:
            #for new_names in new_color_genes:
                #color_genes_handle.write(new_names + "\t" + ",".join([str(ncg) for ncg in new_color_genes[new_names]]) + "\n")
    
    if len(new_color_domains) > 0:
        #if len(new_color_domains) < 10:
            #print("   Saving new color names for domains " + ", ".join(new_color_domains.keys()))
        #else:
            #print("   Saving new color names for 10+ domains")
            
        with open(domains_color_file, "a") as color_domains_handle:
            for new_names in new_color_domains:
                color_domains_handle.write(new_names + "\t" + ",".join([str(ncdom) for ncdom in new_color_domains[new_names]]) + "\n")

    
    mode = "a" if write_html == True else "w"
    with open(outputfile, mode) as handle:
        handle.write(SVG_TEXT)

