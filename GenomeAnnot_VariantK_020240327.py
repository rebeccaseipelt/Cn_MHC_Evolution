#sort quality in vcf file and also exclude parent strain variants

#help on string/integer/float issues with numbers
#https://www.learndatasci.com/solutions/python-valueerror-invalid-literal-int-base-10/

#help on joining with a tab separator
#https://www.geeksforgeeks.org/python-string-concatenation

#help on removing duplicates from a list
#https://www.geeksforgeeks.org/python-ways-to-remove-duplicates-from-list/

#help on converting list to string
#https://www.geeksforgeeks.org/python-program-to-convert-a-list-to-string/

#help on finding intersection of two lists
#https://www.geeksforgeeks.org/python-intersection-two-lists/

#open output file
OUT1 = open("qualityVariants_strain01.txt", "w")
OUT2 = open("qualityVariants_strain02.txt", "w")
OUT3 = open("qualityVariantsParsed_strain01.txt", "w")
OUT4 = open("qualityVariantsParsed_strain02.txt", "w")

#zero out all counts
count_qual_strain01 = 0
count_total_strain01 = 0
count_qual_nonparent01 = 0

count_qual_strain02 = 0
count_total_strain02 = 0
count_qual_nonparent02 = 0

count_qual_parent = 0
count_total_parent = 0

quality = 0
quality2 = 0

#lists
good_qual_list_strain01 = []
good_qual_list_strain02 = []
good_qual_list_parent = []
goodvariant_strain01 = []
goodvariant_strain02 = []

import re

def parent_check(variant):
   common_variant = 0
   for things in good_qual_list_parent:
      if things == variant:
         common_variant +=1
   if common_variant == 0:
      common_variant_result = "no"
   else:
      common_variant_result = "yes"
   return common_variant_result


def intersection_fnct(set1,set2):
   return set(set1).intersection(set2)

def effect(info_column):
   if re.search(r"HIGH", info_column):
      effect_is = "HIGH"
   elif re.search(r"MODERATE", info_column):
      effect_is = "MODERATE"
   elif re.search(r"LOW", info_column):
      effect_is = "LOW"
   else: 
      effect_is = "unknown"
   return effect_is

def list_to_string(stuff):
   #initialize the string
   string_list = " "
   return (string_list.join(stuff))

#open incoming sample files
#open data file, read it in
#use readlines to get each line as a member of a list
FILE_STRAIN01 = open("SNPeff_strain12703_Galaxy.vcf")
FILE_STRAIN02 = open("SNPeff_strain12938_Galaxy.vcf")
FILE_PARENT = open("SNPeff_strainH99W2_Galaxy.vcf")
#FILE_GENOME = open("CnGenesOnly_20230911.csv")

#open output file
#MY_FILE = open("CountGeneFamily_python.txt", "w")

#pull data into lists for processing
list_strain01 = FILE_STRAIN01.readlines()
list_strain02 = FILE_STRAIN02.readlines()
list_parent = FILE_PARENT.readlines()

#loop over parent variants and catch good qual (>=80)

#loop to split family file row and collect 
#columns 1 (number), 2(descriptor)
for step_parent in list_parent:
   if step_parent.startswith("#"):
      header_parent = step_parent
   elif step_parent.startswith("##"):
      count_header+=1
   else:
      count_total_parent+=1
      parent_row = step_parent.split("\t")
      chr_number_parent = parent_row[0]
      location_parent = parent_row[1]
      quality_parent = float(parent_row[5])
      if int(quality_parent) >= 80:
         variant_name_parent = chr_number_parent + "_" + location_parent
         good_qual_list_parent.append(variant_name_parent)
         count_qual_parent+=1
         

#loop to split family file row and collect 
#columns 1 (number), 2(descriptor)
for step_strain01 in list_strain01:
   #count_total_12600+=1
   if step_strain01.startswith("#"):
      header_strain01 = step_strain01
   else:
      count_total_strain01+=1
      row01 = step_strain01.split("\t")
      quality01 = float(row01[5])
      chr_number01 = row01[0]
      location01 = row01[1]
      variant_name01 = chr_number01 + "_" +  location01
      info01 = row01[7]
      #conditional to get good quality
      if int(quality01) >= 80:
         common_variant01 = parent_check(variant_name01)
         count_qual_strain01+=1
         if common_variant01 == "no":
            count_qual_nonparent01+=1
            strainName = "strain01"
            save_this01 = ("\t".join([variant_name01,info01]))
            goodvariant_strain01.append(save_this01)
            #parseInfo(strainName,variant_name01,info01)
            OUT1.write(save_this01 + "\n")

# http://pcingola.github.io/SnpEff/snpeff/inputoutput/#vcf-files
# Info column has a piece of leftover stuff from SNPEFF that shouldn't be in there:
# so splitting at ANN, then at comma (which should break into gene-specific calls)
# each gene-specific call has 15 sections separated by pipe
# section 1 is allele
# section 2 is annotation (effect)
# section 3 is putative (impact)
# section 4 is Gene name

            big_pieces = info01.split("ANN")
            realInfo = big_pieces[1].split(",")
            #print (realInfo)
            for one in realInfo:
                each_piece = one.split("|")
                #print (each_piece[2])
                gene = each_piece[3]
                annotation = each_piece[1]
                impact = each_piece[2]
                save_these = ("\t".join([strainName,variant_name01,gene,annotation,impact]))
                OUT3.write(save_these + "\n")
                gene = ()
                annotation = ()
                impact = ()
            

for step_strain02 in list_strain02:
   #count_total_12600+=1
   if step_strain02.startswith("#"):
      header_strain02 = step_strain02
   else:
      count_total_strain02+=1
      row02 = step_strain02.split("\t")
      quality02 = float(row02[5])
      chr_number02 = row02[0]
      location02 = row02[1]
      variant_name02 = chr_number02 + "_" +  location02
      info02 = row02[7]
      #conditional to get good quality
      if int(quality02) >= 80:
         common_variant02 = parent_check(variant_name02)
         count_qual_strain02+=1
         if common_variant02 == "no":
            count_qual_nonparent02+=1
            strainName2 = "strain02"
            save_this02 = ("\t".join([variant_name02,info02]))
            goodvariant_strain02.append(save_this02)
            OUT2.write(save_this02 + "\n")
            
            big_pieces2 = info02.split("ANN")
            realInfo2 = big_pieces2[1].split(",")
            #print (realInfo)
            for one2 in realInfo2:
                each_piece2 = one2.split("|")
                #print (each_piece[2])
                gene2 = each_piece2[3]
                annotation2 = each_piece2[1]
                impact2 = each_piece2[2]
                save_these2 = ("\t".join([strainName2,variant_name02,gene2,annotation2,impact2]))
                OUT4.write(save_these2 + "\n")
                gene2 = ()
                annotation2 = ()
                impact2 = ()




# regular expression to capture all gene names
#            genes_12602 = re.findall(r"CNAG_.....",info)
#            #remove duplicates from genes_12600
#            no_duplicates_12602 = []
#            [no_duplicates_12602.append(x) for x in genes_12602 if x not in no_duplicates_12602]
            #collapse gene list to a string
#            capture_all_genes2 = list_to_string(no_duplicates_12602)
            
#            save_this = ("\t".join([variant_name,snpEff,capture_all_genes2]))
#            variant_effect_genes_12602.append(save_this)
##            OUT3.write(save_this + "\n")
            
#            for increment in no_duplicates_12602:
#               gene_list_12602.append(increment)
#               save_this = ("\t".join([increment,snpEff,variant_name]))
#               gene_effect_variant_12602.append(save_this)
            #print(increment)

#no_duplicates_final_12602 = []
#[no_duplicates_final_12602.append(x) for x in gene_list_12602 if x not in no_duplicates_final_12602]

#for lines in no_duplicates_final_12602:
#   OUT4.write(lines + "\n")

#common_genes = intersection_fnct(no_duplicates_final_12600, no_duplicates_final_12602)

#for rows in common_genes:
#   OUT5.write(rows + "\n")

#print (no_duplicates_final_12600)

print ("There are "+ str(count_total_parent) + " total variants in parent and " + str(count_qual_parent) + " are good quality." )

print ("\n\nThere are "+ str(count_total_strain01) + " total variants in strain 01 and " + str(count_qual_strain01) + " are good quality.")
print ("There are "+ str(count_qual_nonparent01) + " non-parent good quality variants in strain01")

print ("\n\nThere are "+ str(count_total_strain02) + " total variants in strain02 and " + str(count_qual_strain02) + " are good quality.")
print ("There are "+ str(count_qual_nonparent02) + " non-parent good quality variants in strain02")


#print (*gene_effect_variant_12600, sep = "\n")
#print (*good_qual_list_h, sep = "\n")
#print (*no_duplicates_final_12600, sep = "\n")
#print (*variant_effect_genes_12600, sep = "\n")
