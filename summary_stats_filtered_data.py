import sys

# USE THIS TO GET SUMMARY STATS ON THE AWK FILTERED FILE


# FOR THE individual samples hets and homs ref and homs alt info it prints to screen This file may need to be split on the first column (which is the sample name) and then it can be plot in R.

my_vcf =open(sys.argv[1],'r')
#my_vcf = open("/Users/Maggie/Documents/Duckweed/short5.vcf", "r")


#my_vcf = open("/Users/Maggie/Documents/Duckweed/reduced4.txt", "r")


summ_stats =  open("site_stats_for_filterd_data_dec2.csv", 'w')
#summ_stats =  open("summary_stats_by_site.csv", 'w')
sample_list=['BP1', 'CC1_1', 'CC3_3', 'CC4_1', 'DD7', 'GP10_3', 'GP14_4', 'GP2_3', 'GP4_2', 'GP6_5', 'GP8_1', 'HFA10', 'HFB11', 'ML1_1', 'ML3_1', 'RB1_1', 'RB2_7', 'RB3_7', 'RC1_1', 'RC3_2', 'RC_2_3', 'RD24', 'RL1_1', 'RL2_1', 'RR2_1', 'RR3_6', 'RU100', 'RU102', 'RU103', 'RU186', 'RU195', 'RU206', 'RU408', 'RU410', 'RU448', 'RU99']

def CalcWind(chrom, mid, last_pos, stats_wanted, ind_stats_wanted):  
    summ_stats.write( "\n")
    
    summ_stats.write(chrom +", " + str(mid) + ", " + str(last_pos)+ ", ")
    
    for key, v in stats_wanted.items():
	#print key #keys' are invariant and variant
	# v is the dictionary of stuff for each invar and var
	#print v
	for keys, stuff in v.items():
	    
#	    print stuff #v is the dictionary of coutns for invariant sites and variant sites separately
	    summ_stats.write(str(stuff)+ ", ")



def CalcWind_ind(samp, chrom, mid, ind_stats_wanted):
	#for keys, stuff in v.items():
	    
#	    print stuff #v is the dictionary of coutns for invariant sites and variant sites separately
	 #   summ_stats.write(str(stuff)+ ", ")
    lists= ind_stats_wanted[samp].keys()
    print samp, chrom, mid, ind_stats_wanted[samp][str(lists[0])], ind_stats_wanted[samp][str(lists[1])], ind_stats_wanted[samp][str(lists[2])], ind_stats_wanted[samp][str(lists[3])], ind_stats_wanted[samp][str(lists[4])], ind_stats_wanted[samp][str(lists[5])], ind_stats_wanted[samp][str(lists[6])]
    #print lists
    #I think this only output Het, HOM_R and HOM_ALT calls for filtered sites.
    
    
def count_stats(line, stats, ind_stats_wanted):  
#the first sample info is at sline[5] so I need to make i =2 here    
    i=2
    k=2
    dp=0
    if sline[3] == ".":
        varType = "invariant"
        stats_wanted[varType]['count']+=1
    else:
        varType = "variant"
        stats_wanted[varType]['count']+=1      
	
    samp_num=0
    depth=0
    for sample in sample_list:
        #i is a counter to loop through the genotypes for each sample
        i=i+3
	k=k+3
	#print sline[k]
	if sline[k]!="M" and sline[k+1]!="M":
	    ind_stats_wanted[sample]['Depth']+=int(sline[k+1])
	    ind_stats_wanted[sample]['site_count']+=1
	    
	    ind_stats_wanted[sample]["FILT_"+varType]+=1
	   # print varType,sample, ind_stats_wanted[sample]["FILT_"+varType]
	else: 
	    ind_stats_wanted[sample]['Depth']+=0
#this will count number of invariant sites and number of variant sites
        

        
        if varType == "invariant":
                    #if I type cont it will just go back to start of loop right away
            if (sline[i+1])!="M":
	        samp_num+=1
		depth+=int(sline[i+1]) 
                            
            continue
                #if this is a varaint site then I can get the hets and homs etc.
        else:
            if (sline[i+1])!="M":
	        samp_num+=1
		depth+=int(sline[i+1]) 		
            if sline[i]=="1":
                
                ind_stats_wanted[sample]["Het"]+=1
                stats_wanted[varType]['num_het']+=1     
                #if gentoype is homoz. refernece value is 0
            elif sline[i]=="0":
                ind_stats_wanted[sample]["HOM_R"]+=1                
                stats_wanted[varType]['num_hom_ref']+=1                
	    elif sline[i]=="2":                    #print sample.gt_type
                    #print ind_stats_wanted[sample.sample]
                ind_stats_wanted[sample]["HOM_ALT"]+=1
                stats_wanted[varType]['num_hom_alt']+=1
    stats_wanted[varType]['avg_depth']+=depth/samp_num

#grab all of the individual IDs from the second line in the vcf
numInds = 36
indIDs = ['BP1', 'CC1_1', 'CC3_3', 'CC4_1', 'DD7', 'GP10_3', 'GP14_4', 'GP2_3', 'GP4_2', 'GP6_5', 'GP8_1', 'HFA10', 'HFB11', 'ML1_1', 'ML3_1', 'RB1_1', 'RB2_7', 'RB3_7', 'RC1_1', 'RC3_2', 'RC_2_3', 'RD24', 'RL1_1', 'RL2_1', 'RR2_1', 'RR3_6', 'RU100', 'RU102', 'RU103', 'RU186', 'RU195', 'RU206', 'RU408', 'RU410', 'RU448', 'RU99']


def ResetDicts():
    stats_wanted={'variant':{},'invariant':{}}
    ind_stats_wanted={}
    list_stats=["count", "avg_depth", "no_data", "num_het", "num_hom_alt", "num_hom_ref"]
    
    for x in list_stats:
	#this will make dictionary with keys of the wanted stats
	stats_wanted['variant'][x] = 0.0
	stats_wanted['invariant'][x] = 0.0
    
    ind_stats = ["FILT_variant", "FILT_invariant", "Het", "HOM_R", "HOM_ALT", "Depth", "site_count"]    
    for x in indIDs:
	ind_stats_wanted[x] = {}
	for i in ind_stats:
	    ind_stats_wanted[x][i]=0.0
    
	#print ind_stats_wanted[x]
    return stats_wanted, ind_stats_wanted

stats_wanted, ind_stats_wanted = ResetDicts()
    
window=20000

#I will print column names for the indiv_sample statistics wanted
lists=ind_stats_wanted['BP1'].keys()
#print lists
print "sample_ID", "chrom", "mid_point", str(lists[0]), str(lists[1]), str(lists[2]), str(lists[3]) , str(lists[4]), str(lists[5]), str(lists[6])  

#for samp_names in ind_stats_wanted.items():
#    for keys, i in ind_stats_wanted[samp_names].items():#
#	print i
	




# PREP THe columns for my output file

summ_stats.write("chrom" + ", " + "mid " + ", " + "last_pos" )

#this little bit of code writes out the column names to the summary file
for key, v in stats_wanted.items():
    for keys, stuff in v.items():
	#print keys
	#print stuff #v is the dictionary of coutns for invariant sites and variant sites separately
	summ_stats.write(str(key)+"_"+str(keys)+ ", ")



window=20000
currentChrom="pseudo0"
start=-1
last_pos=0
for line in my_vcf:
    sline=line.split()
    line_POS=int(sline[1])
    #print line_POS, window
    line_CHROM=str(sline[0])
    
    if start==-1:
        start=int(sline[1])
        currentChrom=str(sline[0])
    
    if line_POS >=window or line_CHROM !=currentChrom:
	#if the line pos is within the next window
	
	if line_CHROM !=currentChrom and last_pos<=window and line_POS<(window+20000):
	    #print "not missing windows", sline
	
	    mid=(window-20000)+20000*0.5
	    
	    #note that if there are more positions missing than the window size then the mid point position will not be correct!
	    
	    CalcWind(currentChrom, mid, last_pos, stats_wanted, ind_stats_wanted)
	    window=window+20000

	    for samp in indIDs:
		CalcWind_ind(samp, currentChrom, mid, ind_stats_wanted)	
	
	    stats_wanted, ind_stats_wanted = ResetDicts()
	else:#it means that the position is beyond an additional window
	    #advance and count windows of nothing
	    mid=(window-20000)+20000*0.5
	   
	    CalcWind(currentChrom, mid, last_pos, stats_wanted, ind_stats_wanted)
	    
	    for samp in indIDs:
		CalcWind_ind(samp, currentChrom, mid, ind_stats_wanted)
	   
	    stats_wanted, ind_stats_wanted = ResetDicts()
	    	
	    window=window+20000
	    
	    while line_POS>(window-20000):
		try:
		    mid=(window-20000)+20000*0.5  
		    window=window+20000	  
		    CalcWind(currentChrom, mid,last_pos, stats_wanted, ind_stats_wanted)
	    

		    for samp in indIDs:
			CalcWind_ind(samp, currentChrom, mid, ind_stats_wanted)
			
		except StopIteration:
		    break		    
    count_stats(line, stats_wanted, ind_stats_wanted)	
    last_pos=line_POS
    
#MUST DEAL WITH GETTING last bit of window at the very end of file		
    if line_CHROM !=currentChrom:
	window=20000    
    currentChrom=str(sline[0])
    
      

CalcWind(currentChrom, mid,last_pos, stats_wanted, ind_stats_wanted)

for samp in indIDs:
    CalcWind_ind(samp, currentChrom, mid, ind_stats_wanted)
    
summ_stats.close()

