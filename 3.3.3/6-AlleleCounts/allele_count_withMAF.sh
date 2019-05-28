##### TL - 15/02/16


### This script is an example of a parser to filter SNPs based on the MAF or minimum coverage. 
### It is possible to easily (more easily than below)


## The input file for the script below is a RC file (popoolation) containing 18 populations

#chr   pos     rc      allele_count    allele_states   deletion_sum    snp_type        major_alleles(maa)      minor_alleles(mia)      maa_1   maa_2   maa_3   maa_4   maa_5   maa_6   maa_7   maa_8   mia_1   mia_2   mia_3   mia_4
#   mia_5   mia_6   mia_7   mia_8
#Sc0000000       585     A       2       C/A     0       rc      CNCNNNCC        ANANNNAA        18/30   0/24    20/28   0/17    0/20    0/23    21/36   19/28   10/30   0/24    8/28    0/17    0/20    0/23    13/36   8/28
#Sc0000000       587     A       2       G/A     0       rc      GNGNNNGG        ANANNNAA        17/30   0/24    19/28   0/17    0/20    0/23    24/36   21/28   13/30   0/24    9/28    0/17    0/20    0/23    12/36   7/28


## The output file is a list of filtered variants

# scaffold	pos	ref	Nb_all	all1	all2	sumdeletions	pattern	all1_per_pop	all2_per_pop	counts_all1_pop1	counts_all2_pop1	total_counts_pop1	counts_all1_pop2...
#Sc0000000	603	C	2	C	T	0	pop	CCCCCCCCCCCCCCCCCC	NNNNNNTNNNNNNNTTTN	65	0	65	101	0	101	65	0	65	80      0	80	92	0	92	93	0	93	85	9	94	91	0	91	84	0	84	87	0	87	72	0	72	59	0	59	74	0	74	71	0	71	60	6	66	81	2	83	50	27	77	62	0	62


################ usage ##################
display_usage() { 
	echo -e "This script requires three arguments: " 
	echo -e "\nUSAGE:\n./$0 [path_to_RC_file] [output-prefix] [tmpfile]\n" 
	} 
# check whether user had supplied -h or --help . If yes display usage 
	if [[ ( $1 == "--help") ||  $1 == "-h" ]] 
	then 
		display_usage
		exit 0
	fi 
# if less than two arguments supplied, display usage 
	if [  $# -le 1 ] 
	then 
		echo "WARNING: "
		display_usage
		exit 1
	fi 
#################### parameters
# lire le fichier spécifié par l'utilisateur
### CHANGE TO YOUR PATH###
pwd_inputfile="$1"
tmp00="$3"
prefix="$2"

cutoffMAF=$(echo "0.02") # MAF used for SNPs in the whole dataset


#################### summary of the parameters ###########
startdate=$(date)
output_rappel=$(echo "$prefix""_parameters")
echo "Rappel des paramètres pour simulation du $startdate" > $output_rappel
echo "pwd_infile = $pwd_inputfile <= attention ça doit contenir le path complet au fichier" >> $output_rappel
echo "prefix output files = $prefix" >> $output_rappel
echo "MAF used for filtering = $cutoffMAF" >> $output_rappel

##################### Starting filtering #######################
## step 1: basic transformations ####
#Extract only biallelic sites
pwd_biallelic_sites=$(echo "$prefix""_biallelic_sites")
awk '$4== "2" {print $0}' $pwd_inputfile > $pwd_biallelic_sites

#extract major and minor alleles per locus and per pop + replace "/" by "	"
pwd_biallelic_sites_extended=$(echo "$pwd_biallelic_sites""_extended")
rm $pwd_biallelic_sites_extended
while read line; do echo "$line" > $tmp00; major=$(awk '{print $8}' $tmp00); transfomajor=$(echo "$major" | sed 's/A/A	/g' | sed 's/C/C	/g' | sed 's/G/G	/g' | sed 's/T/T	/g' | sed 's/N/N	/g');  minor=$(awk '{print $9}' $tmp00); transfominor=$(echo "$minor" | sed 's/A/A	/g' | sed 's/C/C	/g' | sed 's/G/G	/g' | sed 's/T/T	/g' | sed 's/N/N	/g');  before=$(echo "$major	$minor"); after=$(echo "$major	$transfomajor$minor	$transfominor"); printable=$(awk '{print $0}' $tmp00 | sed "s/$before/$after/g" | sed 's;/;	;g'); echo "$printable" >> $pwd_biallelic_sites_extended; done < $pwd_biallelic_sites
# step 3: main loops ####
# count reference alleles - attention le nombre de colomnes dépend du nombre de populations, ici uniquement pour 8 pops.
pwd_biallelic_sites_controlbugs=$(echo "$pwd_biallelic_sites""_controlbugs") # contrôle des bugs du script: verif des 8 loops
pwd_biallelic_sites_counts_infosites=$(echo "$pwd_biallelic_sites""_infosites") # info sites générale
pwd_biallelic_sites_counts_infotest=$(echo "$pwd_biallelic_sites""_infotest") # contrôle des bugs du script: verif test fin du script (assez redondant avec infosites, à désactiver par l'utilisateur)
pwd_biallelic_sites_counts=$(echo "$pwd_biallelic_sites""_counts") # sortie standard sous la forme count_ref count_alt totalref+alt
pwd_biallelic_sites_counts2n=$(echo "$pwd_biallelic_sites""_counts2n") #une ligne pour l'allele ref (avec count_ref & totalref+alt), une ligne
pwd_biallelic_sites_counts_excludedMAF=$(echo "$pwd_biallelic_sites""_counts_excludedMAF") # sortie standard sous la forme count_ref count_alt totalref+alt
pwd_biallelic_sites_counts2n_excludedMAF=$(echo "$pwd_biallelic_sites""_counts2n_excludedMAF") #une ligne pour l'allele ref (avec count_ref & totalref+alt), une ligne pour l'allele alt (avec count_alt & totalref+alt) - ça sera probablement plus simple pour faire des analyses sous bayenv par la suite
pwd_biallelic_sites_counts_excludedNA=$(echo "$pwd_biallelic_sites""_counts_excludedNA") # sortie standard sous la forme count_ref count_alt totalref+alt
pwd_biallelic_sites_counts2n_excludedNA=$(echo "$pwd_biallelic_sites""_counts2n_excludedNA") #une ligne pour l'allele ref (avec count_ref & totalref+alt), une ligne pour l'allele alt (avec count_alt & totalref+alt) - ça sera probablement plus simple pour faire des analyses sous bayenv par la suite
rm $pwd_biallelic_sites_controlbugs $pwd_biallelic_sites_counts $pwd_biallelic_sites_counts2n $pwd_biallelic_sites_counts_excludedMAF $pwd_biallelic_sites_counts2n_excludedMAF $pwd_biallelic_sites_counts_infosites $pwd_biallelic_sites_counts_infotest
echo "scaffold	position	allref	allcount	majorall	minorall	sumdeletions	SNPtype	MostFreqAllPerPop	SecondMostFreqAllPerPop	refcount_pop1	altcount_pop1	refaltcount_pop1	refcount_pop2	altcount_pop2	refaltcount_pop2	refcount_pop3	altcount_pop3	refaltcount_pop3	refcount_pop4	altcount_pop4	refaltcount_pop4	refcount_pop5	altcount_po5	refaltcount_pop5	refcount_pop6	altcount_pop6	refaltcount_pop6	refcount_pop7	altcount_pop7	refaltcount_pop7	refcount_pop8	altcount_pop8	refaltcount_pop8	refcount_pop9	altcount_pop9	refaltcount_pop9	refcount_pop10	altcount_pop10	refaltcount_pop10	refcount_pop11	altcount_pop11	refaltcount_pop11	refcount_pop12	altcount_pop12	refaltcount_pop12	refcount_pop13	altcount_pop13	refaltcount_pop13	refcount_pop14	altcount_pop14	refaltcount_pop14	refcount_pop15	altcount_pop15	refaltcount_pop15	refcount_pop16	altcount_pop16	refaltcount_pop16	refcount_pop17	altcount_pop17	refaltcount_pop17	refcount_pop18	altcount_pop18	refaltcount_pop18" > $pwd_biallelic_sites_counts	
while read line; do
	echo "$line" > $tmp00
	ref=$(awk '{print $3}' $tmp00)
	position=$(awk '{print $1"	"$2"	"$3"	"$4"	"$5"	"$6"	"$7"	"$8"	"$9"	"$28}' $tmp00)
	expectedmajor=$(awk '{print $5}' $tmp00)
	expectedminor=$(awk '{print $6}' $tmp00)
	majorall_pop1=$(awk '{print $10}' $tmp00)
	majorall_pop2=$(awk '{print $11}' $tmp00)
	majorall_pop3=$(awk '{print $12}' $tmp00)
	majorall_pop4=$(awk '{print $13}' $tmp00)
	majorall_pop5=$(awk '{print $14}' $tmp00)
	majorall_pop6=$(awk '{print $15}' $tmp00)
	majorall_pop7=$(awk '{print $16}' $tmp00)
	majorall_pop8=$(awk '{print $17}' $tmp00)
	majorall_pop9=$(awk '{print $18}' $tmp00)
	majorall_pop10=$(awk '{print $19}' $tmp00)
	majorall_pop11=$(awk '{print $20}' $tmp00)
	majorall_pop12=$(awk '{print $21}' $tmp00)
	majorall_pop13=$(awk '{print $22}' $tmp00)
	majorall_pop14=$(awk '{print $23}' $tmp00)
	majorall_pop15=$(awk '{print $24}' $tmp00)
	majorall_pop16=$(awk '{print $25}' $tmp00)
	majorall_pop17=$(awk '{print $26}' $tmp00)
	majorall_pop18=$(awk '{print $27}' $tmp00)
	minorall_pop1=$(awk '{print $29}' $tmp00) # décalage d'une colomne pour passer le $minor voir dans le fichier $pwd_biallelic_sites_extended
	minorall_pop2=$(awk '{print $30}' $tmp00)
	minorall_pop3=$(awk '{print $31}' $tmp00)
	minorall_pop4=$(awk '{print $32}' $tmp00)
	minorall_pop5=$(awk '{print $33}' $tmp00)
	minorall_pop6=$(awk '{print $34}' $tmp00)
	minorall_pop7=$(awk '{print $35}' $tmp00)
	minorall_pop8=$(awk '{print $36}' $tmp00)
	minorall_pop9=$(awk '{print $37}' $tmp00)
	minorall_pop10=$(awk '{print $38}' $tmp00)
	minorall_pop11=$(awk '{print $39}' $tmp00)
	minorall_pop12=$(awk '{print $40}' $tmp00)
	minorall_pop13=$(awk '{print $41}' $tmp00)
	minorall_pop14=$(awk '{print $42}' $tmp00)
	minorall_pop15=$(awk '{print $43}' $tmp00)
	minorall_pop16=$(awk '{print $44}' $tmp00)
	minorall_pop17=$(awk '{print $45}' $tmp00)
	minorall_pop18=$(awk '{print $46}' $tmp00)
	#echo "$ref	$expectedminor	$expectedmajor	$majorall_pop1	$minorall_pop1	$majorall_pop2	$minorpop2 $majorall_pop3... $position"
	count_ref_pop1=$(echo "")
	count_alt_pop1=$(echo "")
	count_ref_alt_pop1=$(echo "")
	count_ref_pop2=$(echo "")
	count_alt_pop2=$(echo "")
	count_ref_alt_pop2=$(echo "")
	count_ref_pop3=$(echo "")
	count_alt_pop3=$(echo "")
	count_ref_alt_pop3=$(echo "")
	count_ref_pop4=$(echo "")
	count_alt_pop4=$(echo "")
	count_ref_alt_pop4=$(echo "")
	count_ref_pop5=$(echo "")
	count_alt_pop5=$(echo "")
	count_ref_alt_pop5=$(echo "")
	count_ref_pop6=$(echo "")
	count_alt_pop6=$(echo "")
	count_ref_alt_pop6=$(echo "")
	count_ref_pop7=$(echo "")
	count_alt_pop7=$(echo "")
	count_ref_alt_pop7=$(echo "")
	count_ref_pop8=$(echo "")
	count_alt_pop8=$(echo "")
	count_ref_alt_pop8=$(echo "")
	count_ref_pop9=$(echo "")
	count_alt_pop9=$(echo "")
	count_ref_alt_pop9=$(echo "")
	count_ref_pop10=$(echo "")
	count_alt_pop10=$(echo "")
	count_ref_alt_pop10=$(echo "")
	count_ref_pop11=$(echo "")
	count_alt_pop11=$(echo "")
	count_ref_alt_pop11=$(echo "")
	count_ref_pop12=$(echo "")
	count_alt_pop12=$(echo "")
	count_ref_alt_pop12=$(echo "")
	count_ref_pop13=$(echo "")
	count_alt_pop13=$(echo "")
	count_ref_alt_pop13=$(echo "")
	count_ref_pop14=$(echo "")
	count_alt_pop14=$(echo "")
	count_ref_alt_pop14=$(echo "")
	count_ref_pop15=$(echo "")
	count_alt_pop15=$(echo "")
	count_ref_alt_pop15=$(echo "")
	count_ref_pop16=$(echo "")
	count_alt_pop16=$(echo "")
	count_ref_alt_pop16=$(echo "")
	count_ref_pop17=$(echo "")
	count_alt_pop17=$(echo "")
	count_ref_alt_pop17=$(echo "")
	count_ref_pop18=$(echo "")
	count_alt_pop18=$(echo "")
	count_ref_alt_pop18=$(echo "")
	# loop pop1
	if [ "$minorall_pop1" == "$expectedminor" ] || [ "$minorall_pop1" == "$expectedmajor" ] || [ "$minorall_pop1" == "N" ]; then
		if [ "$ref" == "$majorall_pop1" ]; then
			count_ref_pop1=$(awk '{print $47}' $tmp00)
			count_alt_pop1=$(awk '{print $83}' $tmp00)
			count_ref_alt_pop1=$(echo "$count_ref_pop1 + $count_alt_pop1" | bc)
		elif [ "$ref" == "$minorall_pop1" ]; then
			count_ref_pop1=$(awk '{print $83}' $tmp00)
			count_alt_pop1=$(awk '{print $47}' $tmp00)
			count_ref_alt_pop1=$(echo "$count_ref_pop1 + $count_alt_pop1" | bc)
		elif [ "$majorall_pop1" == "$minorall_pop1" ]; then
			count_ref_pop1=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop1=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop1=$(echo "NA") 
			echo "cas 2: majorall_pop1 == minorall_pop1 pour $majorall_pop1 == $minorall_pop1 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop1" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop1 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop1" ]; then
				count_ref_pop1=$(awk '{print $47}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop1=$(awk '{print $83}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop1=$(echo "$count_ref_pop1 + $count_alt_pop1" | bc)
			else
				count_ref_pop1=$(awk '{print $83}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop1=$(awk '{print $47}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop1=$(echo "$count_ref_pop1 + $count_alt_pop1" | bc)
			fi
		else
			count_ref_pop1=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop1=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop1=$(echo "NA")
			echo "autre cas à analyser: majorall_pop1=$majorall_pop1 minorall_pop1=$minorall_pop1 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop1=$(echo "NA") 
		count_alt_pop1=$(echo "NA")
		count_ref_alt_pop1=$(echo "NA")
		echo "cas 1: minorall_pop1 == $minorall_pop1 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	# loop pop2
	if [ "$minorall_pop2" == "$expectedminor" ] || [ "$minorall_pop2" == "$expectedmajor" ] || [ "$minorall_pop2" == "N" ]; then
		if [ "$ref" == "$majorall_pop2" ]; then
			count_ref_pop2=$(awk '{print $49}' $tmp00)
			count_alt_pop2=$(awk '{print $85}' $tmp00)
			count_ref_alt_pop2=$(echo "$count_ref_pop2 + $count_alt_pop2" | bc)
		elif [ "$ref" == "$minorall_pop2" ]; then
			count_ref_pop2=$(awk '{print $85}' $tmp00)
			count_alt_pop2=$(awk '{print $49}' $tmp00)
			count_ref_alt_pop2=$(echo "$count_ref_pop2 + $count_alt_pop2" | bc)
		elif [ "$majorall_pop2" == "$minorall_pop2" ]; then
			count_ref_pop2=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop2=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop2=$(echo "NA") 
			echo "cas 2: majorall_pop2 == minorall_pop2 pour $majorall_pop2 == $minorall_pop2 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop2" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop2 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop2" ]; then
				count_ref_pop2=$(awk '{print $49}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop2=$(awk '{print $85}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop2=$(echo "$count_ref_pop2 + $count_alt_pop2" | bc)
			else
				count_ref_pop2=$(awk '{print $85}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop2=$(awk '{print $49}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop2=$(echo "$count_ref_pop2 + $count_alt_pop2" | bc)
			fi
		else
			count_ref_pop2=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop2=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop2=$(echo "NA")
			echo "autre cas à analyser: majorall_pop2=$majorall_pop2 minorall_pop2=$minorall_pop2 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop2=$(echo "NA") 
		count_alt_pop2=$(echo "NA")
		count_ref_alt_pop2=$(echo "NA")
		echo "cas 1: minorall_pop2 == $minorall_pop2 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	# loop pop3
	if [ "$minorall_pop3" == "$expectedminor" ] || [ "$minorall_pop3" == "$expectedmajor" ] || [ "$minorall_pop3" == "N" ]; then
		if [ "$ref" == "$majorall_pop3" ]; then
			count_ref_pop3=$(awk '{print $51}' $tmp00)
			count_alt_pop3=$(awk '{print $87}' $tmp00)
			count_ref_alt_pop3=$(echo "$count_ref_pop3 + $count_alt_pop3" | bc)
		elif [ "$ref" == "$minorall_pop3" ]; then
			count_ref_pop3=$(awk '{print $87}' $tmp00)
			count_alt_pop3=$(awk '{print $51}' $tmp00)
			count_ref_alt_pop3=$(echo "$count_ref_pop3 + $count_alt_pop3" | bc)
		elif [ "$majorall_pop3" == "$minorall_pop3" ]; then
			count_ref_pop3=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop3=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop3=$(echo "NA") 
			echo "cas 2: majorall_pop3 == minorall_pop3 pour $majorall_pop3 == $minorall_pop3 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop3" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop3 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop3" ]; then
				count_ref_pop3=$(awk '{print $51}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop3=$(awk '{print $87}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop3=$(echo "$count_ref_pop3 + $count_alt_pop3" | bc)
			else
				count_ref_pop3=$(awk '{print $87}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop3=$(awk '{print $51}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop3=$(echo "$count_ref_pop3 + $count_alt_pop3" | bc)
			fi
		else
			count_ref_pop3=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop3=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop3=$(echo "NA")
			echo "autre cas à analyser: majorall_pop3=$majorall_pop3 minorall_pop3=$minorall_pop3 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop3=$(echo "NA") 
		count_alt_pop3=$(echo "NA")
		count_ref_alt_pop3=$(echo "NA")
		echo "cas 1: minorall_pop3 == $minorall_pop3 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	# loop pop4
	if [ "$minorall_pop4" == "$expectedminor" ] || [ "$minorall_pop4" == "$expectedmajor" ] || [ "$minorall_pop4" == "N" ]; then
		if [ "$ref" == "$majorall_pop4" ]; then
			count_ref_pop4=$(awk '{print $53}' $tmp00)
			count_alt_pop4=$(awk '{print $89}' $tmp00)
			count_ref_alt_pop4=$(echo "$count_ref_pop4 + $count_alt_pop4" | bc)
		elif [ "$ref" == "$minorall_pop4" ]; then
			count_ref_pop4=$(awk '{print $89}' $tmp00)
			count_alt_pop4=$(awk '{print $53}' $tmp00)
			count_ref_alt_pop4=$(echo "$count_ref_pop4 + $count_alt_pop4" | bc)
		elif [ "$majorall_pop4" == "$minorall_pop4" ]; then
			count_ref_pop4=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop4=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop4=$(echo "NA") 
			echo "cas 2: majorall_pop4 == minorall_pop4 pour $majorall_pop4 == $minorall_pop4 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop4" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop4 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop4" ]; then
				count_ref_pop4=$(awk '{print $53}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop4=$(awk '{print $89}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop4=$(echo "$count_ref_pop4 + $count_alt_pop4" | bc)
			else
				count_ref_pop4=$(awk '{print $89}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop4=$(awk '{print $53}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop4=$(echo "$count_ref_pop4 + $count_alt_pop4" | bc)
			fi
		else
			count_ref_pop4=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop4=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop4=$(echo "NA")
			echo "autre cas à analyser: majorall_pop4=$majorall_pop4 minorall_pop4=$minorall_pop4 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop4=$(echo "NA") 
		count_alt_pop4=$(echo "NA")
		count_ref_alt_pop4=$(echo "NA")
		echo "cas 1: minorall_pop4 == $minorall_pop4 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	# loop pop5
	if [ "$minorall_pop5" == "$expectedminor" ] || [ "$minorall_pop5" == "$expectedmajor" ] || [ "$minorall_pop5" == "N" ]; then
		if [ "$ref" == "$majorall_pop5" ]; then
			count_ref_pop5=$(awk '{print $55}' $tmp00)
			count_alt_pop5=$(awk '{print $91}' $tmp00)
			count_ref_alt_pop5=$(echo "$count_ref_pop5 + $count_alt_pop5" | bc)
		elif [ "$ref" == "$minorall_pop5" ]; then
			count_ref_pop5=$(awk '{print $91}' $tmp00)
			count_alt_pop5=$(awk '{print $55}' $tmp00)
			count_ref_alt_pop5=$(echo "$count_ref_pop5 + $count_alt_pop5" | bc)
		elif [ "$majorall_pop5" == "$minorall_pop5" ]; then
			count_ref_pop5=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop5=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop5=$(echo "NA") 
			echo "cas 2: majorall_pop5 == minorall_pop5 pour $majorall_pop5 == $minorall_pop5 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop5" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop5 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop5" ]; then
				count_ref_pop5=$(awk '{print $55}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop5=$(awk '{print $91}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop5=$(echo "$count_ref_pop5 + $count_alt_pop5" | bc)
			else
				count_ref_pop5=$(awk '{print $91}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop5=$(awk '{print $55}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop5=$(echo "$count_ref_pop5 + $count_alt_pop5" | bc)
			fi
		else
			count_ref_pop5=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop5=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop5=$(echo "NA")
			echo "autre cas à analyser: majorall_pop5=$majorall_pop5 minorall_pop5=$minorall_pop5 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop5=$(echo "NA") 
		count_alt_pop5=$(echo "NA")
		count_ref_alt_pop5=$(echo "NA")
		echo "cas 1: minorall_pop5 == $minorall_pop5 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	# loop pop6
	if [ "$minorall_pop6" == "$expectedminor" ] || [ "$minorall_pop6" == "$expectedmajor" ] || [ "$minorall_pop6" == "N" ]; then
		if [ "$ref" == "$majorall_pop6" ]; then
			count_ref_pop6=$(awk '{print $57}' $tmp00)
			count_alt_pop6=$(awk '{print $93}' $tmp00)
			count_ref_alt_pop6=$(echo "$count_ref_pop6 + $count_alt_pop6" | bc)
		elif [ "$ref" == "$minorall_pop6" ]; then
			count_ref_pop6=$(awk '{print $93}' $tmp00)
			count_alt_pop6=$(awk '{print $57}' $tmp00)
			count_ref_alt_pop6=$(echo "$count_ref_pop6 + $count_alt_pop6" | bc)
		elif [ "$majorall_pop6" == "$minorall_pop6" ]; then
			count_ref_pop6=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop6=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop6=$(echo "NA") 
			echo "cas 2: majorall_pop6 == minorall_pop6 pour $majorall_pop6 == $minorall_pop6 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop6" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop6 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop6" ]; then
				count_ref_pop6=$(awk '{print $57}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop6=$(awk '{print $93}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop6=$(echo "$count_ref_pop6 + $count_alt_pop6" | bc)
			else
				count_ref_pop6=$(awk '{print $93}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop6=$(awk '{print $57}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop6=$(echo "$count_ref_pop6 + $count_alt_pop6" | bc)
			fi
		else
			count_ref_pop6=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop6=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop6=$(echo "NA")
			echo "autre cas à analyser: majorall_pop6=$majorall_pop6 minorall_pop6=$minorall_pop6 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop6=$(echo "NA") 
		count_alt_pop6=$(echo "NA")
		count_ref_alt_pop6=$(echo "NA")
		echo "cas 1: minorall_pop6 == $minorall_pop6 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop7
	if [ "$minorall_pop7" == "$expectedminor" ] || [ "$minorall_pop7" == "$expectedmajor" ] || [ "$minorall_pop7" == "N" ]; then
		if [ "$ref" == "$majorall_pop7" ]; then
			count_ref_pop7=$(awk '{print $59}' $tmp00)
			count_alt_pop7=$(awk '{print $95}' $tmp00)
			count_ref_alt_pop7=$(echo "$count_ref_pop7 + $count_alt_pop7" | bc)
		elif [ "$ref" == "$minorall_pop7" ]; then
			count_ref_pop7=$(awk '{print $95}' $tmp00)
			count_alt_pop7=$(awk '{print $59}' $tmp00)
			count_ref_alt_pop7=$(echo "$count_ref_pop7 + $count_alt_pop7" | bc)
		elif [ "$majorall_pop7" == "$minorall_pop7" ]; then
			count_ref_pop7=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop7=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop7=$(echo "NA") 
			echo "cas 2: majorall_pop7 == minorall_pop7 pour $majorall_pop7 == $minorall_pop7 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop7" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop7 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop7" ]; then
				count_ref_pop7=$(awk '{print $59}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop7=$(awk '{print $95}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop7=$(echo "$count_ref_pop7 + $count_alt_pop7" | bc)
			else
				count_ref_pop7=$(awk '{print $95}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop7=$(awk '{print $59}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop7=$(echo "$count_ref_pop7 + $count_alt_pop7" | bc)
			fi
		else
			count_ref_pop7=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop7=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop7=$(echo "NA")
			echo "autre cas à analyser: majorall_pop7=$majorall_pop7 minorall_pop7=$minorall_pop7 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop7=$(echo "NA") 
		count_alt_pop7=$(echo "NA")
		count_ref_alt_pop7=$(echo "NA")
		echo "cas 1: minorall_pop7 == $minorall_pop7 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop8
	if [ "$minorall_pop8" == "$expectedminor" ] || [ "$minorall_pop8" == "$expectedmajor" ] || [ "$minorall_pop8" == "N" ]; then
		if [ "$ref" == "$majorall_pop8" ]; then
			count_ref_pop8=$(awk '{print $61}' $tmp00)
			count_alt_pop8=$(awk '{print $97}' $tmp00)
			count_ref_alt_pop8=$(echo "$count_ref_pop8 + $count_alt_pop8" | bc)
		elif [ "$ref" == "$minorall_pop8" ]; then
			count_ref_pop8=$(awk '{print $97}' $tmp00)
			count_alt_pop8=$(awk '{print $61}' $tmp00)
			count_ref_alt_pop8=$(echo "$count_ref_pop8 + $count_alt_pop8" | bc)
		elif [ "$majorall_pop8" == "$minorall_pop8" ]; then
			count_ref_pop8=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop8=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop8=$(echo "NA") 
			echo "cas 2: majorall_pop8 == minorall_pop8 pour $majorall_pop8 == $minorall_pop8 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop8" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop8 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop8" ]; then
				count_ref_pop8=$(awk '{print $61}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop8=$(awk '{print $97}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop8=$(echo "$count_ref_pop8 + $count_alt_pop8" | bc)
			else
				count_ref_pop8=$(awk '{print $97}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop8=$(awk '{print $61}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop8=$(echo "$count_ref_pop8 + $count_alt_pop8" | bc)
			fi
		else
			count_ref_pop8=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop8=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop8=$(echo "NA")
			echo "autre cas à analyser: majorall_pop8=$majorall_pop8 minorall_pop8=$minorall_pop8 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop8=$(echo "NA") 
		count_alt_pop8=$(echo "NA")
		count_ref_alt_pop8=$(echo "NA")
		echo "cas 1: minorall_pop8 == $minorall_pop8 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop9
	if [ "$minorall_pop9" == "$expectedminor" ] || [ "$minorall_pop9" == "$expectedmajor" ] || [ "$minorall_pop9" == "N" ]; then
		if [ "$ref" == "$majorall_pop9" ]; then
			count_ref_pop9=$(awk '{print $63}' $tmp00)
			count_alt_pop9=$(awk '{print $99}' $tmp00)
			count_ref_alt_pop9=$(echo "$count_ref_pop9 + $count_alt_pop9" | bc)
		elif [ "$ref" == "$minorall_pop9" ]; then
			count_ref_pop9=$(awk '{print $99}' $tmp00)
			count_alt_pop9=$(awk '{print $63}' $tmp00)
			count_ref_alt_pop9=$(echo "$count_ref_pop9 + $count_alt_pop9" | bc)
		elif [ "$majorall_pop9" == "$minorall_pop9" ]; then
			count_ref_pop9=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop9=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop9=$(echo "NA") 
			echo "cas 2: majorall_pop9 == minorall_pop9 pour $majorall_pop9 == $minorall_pop9 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop9" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop9 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop9" ]; then
				count_ref_pop9=$(awk '{print $63}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop9=$(awk '{print $99}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop9=$(echo "$count_ref_pop9 + $count_alt_pop9" | bc)
			else
				count_ref_pop9=$(awk '{print $99}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop9=$(awk '{print $63}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop9=$(echo "$count_ref_pop9 + $count_alt_pop9" | bc)
			fi
		else
			count_ref_pop9=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop9=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop9=$(echo "NA")
			echo "autre cas à analyser: majorall_pop9=$majorall_pop9 minorall_pop9=$minorall_pop9 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop9=$(echo "NA") 
		count_alt_pop9=$(echo "NA")
		count_ref_alt_pop9=$(echo "NA")
		echo "cas 1: minorall_pop9 == $minorall_pop9 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop10
	if [ "$minorall_pop10" == "$expectedminor" ] || [ "$minorall_pop10" == "$expectedmajor" ] || [ "$minorall_pop10" == "N" ]; then
		if [ "$ref" == "$majorall_pop10" ]; then
			count_ref_pop10=$(awk '{print $65}' $tmp00)
			count_alt_pop10=$(awk '{print $101}' $tmp00)
			count_ref_alt_pop10=$(echo "$count_ref_pop10 + $count_alt_pop10" | bc)
		elif [ "$ref" == "$minorall_pop10" ]; then
			count_ref_pop10=$(awk '{print $101}' $tmp00)
			count_alt_pop10=$(awk '{print $65}' $tmp00)
			count_ref_alt_pop10=$(echo "$count_ref_pop10 + $count_alt_pop10" | bc)
		elif [ "$majorall_pop10" == "$minorall_pop10" ]; then
			count_ref_pop10=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop10=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop10=$(echo "NA") 
			echo "cas 2: majorall_pop10 == minorall_pop10 pour $majorall_pop10 == $minorall_pop10 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop10" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop10 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop10" ]; then
				count_ref_pop10=$(awk '{print $65}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop10=$(awk '{print $101}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop10=$(echo "$count_ref_pop10 + $count_alt_pop10" | bc)
			else
				count_ref_pop10=$(awk '{print $101}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop10=$(awk '{print $65}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop10=$(echo "$count_ref_pop10 + $count_alt_pop10" | bc)
			fi
		else
			count_ref_pop10=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop10=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop10=$(echo "NA")
			echo "autre cas à analyser: majorall_pop10=$majorall_pop10 minorall_pop10=$minorall_pop10 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop10=$(echo "NA") 
		count_alt_pop10=$(echo "NA")
		count_ref_alt_pop10=$(echo "NA")
		echo "cas 1: minorall_pop10 == $minorall_pop10 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop11
	if [ "$minorall_pop11" == "$expectedminor" ] || [ "$minorall_pop11" == "$expectedmajor" ] || [ "$minorall_pop11" == "N" ]; then
		if [ "$ref" == "$majorall_pop11" ]; then
			count_ref_pop11=$(awk '{print $67}' $tmp00)
			count_alt_pop11=$(awk '{print $103}' $tmp00)
			count_ref_alt_pop11=$(echo "$count_ref_pop11 + $count_alt_pop11" | bc)
		elif [ "$ref" == "$minorall_pop11" ]; then
			count_ref_pop11=$(awk '{print $103}' $tmp00)
			count_alt_pop11=$(awk '{print $67}' $tmp00)
			count_ref_alt_pop11=$(echo "$count_ref_pop11 + $count_alt_pop11" | bc)
		elif [ "$majorall_pop11" == "$minorall_pop11" ]; then
			count_ref_pop11=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop11=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop11=$(echo "NA") 
			echo "cas 2: majorall_pop11 == minorall_pop11 pour $majorall_pop11 == $minorall_pop11 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop11" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop11 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop11" ]; then
				count_ref_pop11=$(awk '{print $67}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop11=$(awk '{print $103}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop11=$(echo "$count_ref_pop11 + $count_alt_pop11" | bc)
			else
				count_ref_pop11=$(awk '{print $103}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop11=$(awk '{print $67}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop11=$(echo "$count_ref_pop11 + $count_alt_pop11" | bc)
			fi
		else
			count_ref_pop11=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop11=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop1=$(echo "NA")
			echo "autre cas à analyser: majorall_pop7=$majorall_pop11 minorall_pop7=$minorall_pop11 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop11=$(echo "NA") 
		count_alt_pop11=$(echo "NA")
		count_ref_alt_pop11=$(echo "NA")
		echo "cas 1: minorall_pop11 == $minorall_pop11 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop12
	if [ "$minorall_pop12" == "$expectedminor" ] || [ "$minorall_pop12" == "$expectedmajor" ] || [ "$minorall_pop12" == "N" ]; then
		if [ "$ref" == "$majorall_pop12" ]; then
			count_ref_pop12=$(awk '{print $69}' $tmp00)
			count_alt_pop12=$(awk '{print $105}' $tmp00)
			count_ref_alt_pop12=$(echo "$count_ref_pop12 + $count_alt_pop12" | bc)
		elif [ "$ref" == "$minorall_pop12" ]; then
			count_ref_pop12=$(awk '{print $105}' $tmp00)
			count_alt_pop12=$(awk '{print $69}' $tmp00)
			count_ref_alt_pop12=$(echo "$count_ref_pop12 + $count_alt_pop12" | bc)
		elif [ "$majorall_pop12" == "$minorall_pop12" ]; then
			count_ref_pop12=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop12=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop12=$(echo "NA") 
			echo "cas 2: majorall_pop12 == minorall_pop12 pour $majorall_pop12 == $minorall_pop12 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop12" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop12 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop12" ]; then
				count_ref_pop12=$(awk '{print $69}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop12=$(awk '{print $105}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop12=$(echo "$count_ref_pop12 + $count_alt_pop12" | bc)
			else
				count_ref_pop12=$(awk '{print $105}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop12=$(awk '{print $69}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop12=$(echo "$count_ref_pop12 + $count_alt_pop12" | bc)
			fi
		else
			count_ref_pop12=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop12=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop12=$(echo "NA")
			echo "autre cas à analyser: majorall_pop12=$majorall_pop12 minorall_pop12=$minorall_pop12 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop12=$(echo "NA") 
		count_alt_pop12=$(echo "NA")
		count_ref_alt_pop12=$(echo "NA")
		echo "cas 1: minorall_pop12 == $minorall_pop12 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop13
	if [ "$minorall_pop13" == "$expectedminor" ] || [ "$minorall_pop13" == "$expectedmajor" ] || [ "$minorall_pop13" == "N" ]; then
		if [ "$ref" == "$majorall_pop13" ]; then
			count_ref_pop13=$(awk '{print $71}' $tmp00)
			count_alt_pop13=$(awk '{print $107}' $tmp00)
			count_ref_alt_pop13=$(echo "$count_ref_pop13 + $count_alt_pop13" | bc)
		elif [ "$ref" == "$minorall_pop13" ]; then
			count_ref_pop13=$(awk '{print $107}' $tmp00)
			count_alt_pop13=$(awk '{print $71}' $tmp00)
			count_ref_alt_pop13=$(echo "$count_ref_pop13 + $count_alt_pop13" | bc)
		elif [ "$majorall_pop13" == "$minorall_pop13" ]; then
			count_ref_pop13=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop13=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop13=$(echo "NA") 
			echo "cas 2: majorall_pop13 == minorall_pop13 pour $majorall_pop13 == $minorall_pop13 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop13" == "N" ]; then
			echo "cas 3: minorall_pop13 == N où $majorall_pop13 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop13" ]; then
				count_ref_pop13=$(awk '{print $71}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop13=$(awk '{print $107}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop13=$(echo "$count_ref_pop13 + $count_alt_pop13" | bc)
			else
				count_ref_pop13=$(awk '{print $107}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop13=$(awk '{print $71}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop13=$(echo "$count_ref_pop13 + $count_alt_pop13" | bc)
			fi
		else
			count_ref_pop13=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop13=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop13=$(echo "NA")
			echo "autre cas à analyser: majorall_pop13=$majorall_pop13 minorall_pop13=$minorall_pop13 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop13=$(echo "NA") 
		count_alt_pop13=$(echo "NA")
		count_ref_alt_pop13=$(echo "NA")
		echo "cas 1: minorall_pop13 == $minorall_pop13 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop14
	if [ "$minorall_pop14" == "$expectedminor" ] || [ "$minorall_pop14" == "$expectedmajor" ] || [ "$minorall_pop14" == "N" ]; then
		if [ "$ref" == "$majorall_pop14" ]; then
			count_ref_pop14=$(awk '{print $73}' $tmp00)
			count_alt_pop14=$(awk '{print $109}' $tmp00)
			count_ref_alt_pop14=$(echo "$count_ref_pop14 + $count_alt_pop14" | bc)
		elif [ "$ref" == "$minorall_pop14" ]; then
			count_ref_pop14=$(awk '{print $109}' $tmp00)
			count_alt_pop14=$(awk '{print $73}' $tmp00)
			count_ref_alt_pop14=$(echo "$count_ref_pop14 + $count_alt_pop14" | bc)
		elif [ "$majorall_pop14" == "$minorall_pop14" ]; then
			count_ref_pop14=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop14=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop14=$(echo "NA") 
			echo "cas 2: majorall_pop14 == minorall_pop14 pour $majorall_pop14 == $minorall_pop14 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop14" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop14 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop14" ]; then
				count_ref_pop14=$(awk '{print $73}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop14=$(awk '{print $109}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop14=$(echo "$count_ref_pop14 + $count_alt_pop14" | bc)
			else
				count_ref_pop14=$(awk '{print $109}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop14=$(awk '{print $73}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop14=$(echo "$count_ref_pop14 + $count_alt_pop14" | bc)
			fi
		else
			count_ref_pop14=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop14=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop14=$(echo "NA")
			echo "autre cas à analyser: majorall_pop14=$majorall_pop14 minorall_pop14=$minorall_pop14 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop14=$(echo "NA") 
		count_alt_pop14=$(echo "NA")
		count_ref_alt_pop14=$(echo "NA")
		echo "cas 1: minorall_pop14 == $minorall_pop14 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop15
	if [ "$minorall_pop15" == "$expectedminor" ] || [ "$minorall_pop15" == "$expectedmajor" ] || [ "$minorall_pop15" == "N" ]; then
		if [ "$ref" == "$majorall_pop15" ]; then
			count_ref_pop15=$(awk '{print $75}' $tmp00)
			count_alt_pop15=$(awk '{print $111}' $tmp00)
			count_ref_alt_pop15=$(echo "$count_ref_pop15 + $count_alt_pop15" | bc)
		elif [ "$ref" == "$minorall_pop15" ]; then
			count_ref_pop15=$(awk '{print $111}' $tmp00)
			count_alt_pop15=$(awk '{print $75}' $tmp00)
			count_ref_alt_pop15=$(echo "$count_ref_pop15 + $count_alt_pop15" | bc)
		elif [ "$majorall_pop15" == "$minorall_pop15" ]; then
			count_ref_pop15=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop15=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop15=$(echo "NA") 
			echo "cas 2: majorall_pop15 == minorall_pop15 pour $majorall_pop15 == $minorall_pop15 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop15" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop15 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop15" ]; then
				count_ref_pop15=$(awk '{print $75}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop15=$(awk '{print $111}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop15=$(echo "$count_ref_pop15 + $count_alt_pop15" | bc)
			else
				count_ref_pop15=$(awk '{print $111}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop15=$(awk '{print $75}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop15=$(echo "$count_ref_pop15 + $count_alt_pop15" | bc)
			fi
		else
			count_ref_pop15=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop15=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop15=$(echo "NA")
			echo "autre cas à analyser: majorall_pop15=$majorall_pop15 minorall_pop15=$minorall_pop15 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop15=$(echo "NA") 
		count_alt_pop15=$(echo "NA")
		count_ref_alt_pop15=$(echo "NA")
		echo "cas 1: minorall_pop15 == $minorall_pop15 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop16
	if [ "$minorall_pop16" == "$expectedminor" ] || [ "$minorall_pop16" == "$expectedmajor" ] || [ "$minorall_pop16" == "N" ]; then
		if [ "$ref" == "$majorall_pop16" ]; then
			count_ref_pop16=$(awk '{print $77}' $tmp00)
			count_alt_pop16=$(awk '{print $113}' $tmp00)
			count_ref_alt_pop16=$(echo "$count_ref_pop16 + $count_alt_pop16" | bc)
		elif [ "$ref" == "$minorall_pop16" ]; then
			count_ref_pop16=$(awk '{print $113}' $tmp00)
			count_alt_pop16=$(awk '{print $77}' $tmp00)
			count_ref_alt_pop16=$(echo "$count_ref_pop16 + $count_alt_pop16" | bc)
		elif [ "$majorall_pop16" == "$minorall_pop16" ]; then
			count_ref_pop16=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop16=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop16=$(echo "NA") 
			echo "cas 2: majorall_pop16 == minorall_pop16 pour $majorall_pop16 == $minorall_pop16 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop16" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop16 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop16" ]; then
				count_ref_pop16=$(awk '{print $77}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop16=$(awk '{print $113}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop16=$(echo "$count_ref_pop16 + $count_alt_pop16" | bc)
			else
				count_ref_pop16=$(awk '{print $113}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop16=$(awk '{print $77}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop16=$(echo "$count_ref_pop16 + $count_alt_pop16" | bc)
			fi
		else
			count_ref_pop16=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop16=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop16=$(echo "NA")
			echo "autre cas à analyser: majorall_pop16=$majorall_pop16 minorall_pop16=$minorall_pop16 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop16=$(echo "NA") 
		count_alt_pop16=$(echo "NA")
		count_ref_alt_pop16=$(echo "NA")
		echo "cas 1: minorall_pop16 == $minorall_pop16 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop17
	if [ "$minorall_pop17" == "$expectedminor" ] || [ "$minorall_pop17" == "$expectedmajor" ] || [ "$minorall_pop17" == "N" ]; then
		if [ "$ref" == "$majorall_pop17" ]; then
			count_ref_pop17=$(awk '{print $79}' $tmp00)
			count_alt_pop17=$(awk '{print $115}' $tmp00)
			count_ref_alt_pop17=$(echo "$count_ref_pop17 + $count_alt_pop17" | bc)
		elif [ "$ref" == "$minorall_pop17" ]; then
			count_ref_pop17=$(awk '{print $115}' $tmp00)
			count_alt_pop17=$(awk '{print $79}' $tmp00)
			count_ref_alt_pop17=$(echo "$count_ref_pop17 + $count_alt_pop17" | bc)
		elif [ "$majorall_pop17" == "$minorall_pop17" ]; then
			count_ref_pop17=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop17=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop17=$(echo "NA") 
			echo "cas 2: majorall_pop17 == minorall_pop17 pour $majorall_pop17 == $minorall_pop17 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop17" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop17 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop17" ]; then
				count_ref_pop17=$(awk '{print $79}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop17=$(awk '{print $115}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop17=$(echo "$count_ref_pop17 + $count_alt_pop17" | bc)
			else
				count_ref_pop17=$(awk '{print $115}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop17=$(awk '{print $79}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop17=$(echo "$count_ref_pop17 + $count_alt_pop17" | bc)
			fi
		else
			count_ref_pop17=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop17=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop17=$(echo "NA")
			echo "autre cas à analyser: majorall_pop17=$majorall_pop17 minorall_pop17=$minorall_pop17 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop17=$(echo "NA") 
		count_alt_pop17=$(echo "NA")
		count_ref_alt_pop17=$(echo "NA")
		echo "cas 1: minorall_pop17 == $minorall_pop17 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	#loop pop18
	if [ "$minorall_pop18" == "$expectedminor" ] || [ "$minorall_pop18" == "$expectedmajor" ] || [ "$minorall_pop18" == "N" ]; then
		if [ "$ref" == "$majorall_pop18" ]; then
			count_ref_pop18=$(awk '{print $81}' $tmp00)
			count_alt_pop18=$(awk '{print $117}' $tmp00)
			count_ref_alt_pop18=$(echo "$count_ref_pop18 + $count_alt_pop18" | bc)
		elif [ "$ref" == "$minorall_pop18" ]; then
			count_ref_pop18=$(awk '{print $117}' $tmp00)
			count_alt_pop18=$(awk '{print $81}' $tmp00)
			count_ref_alt_pop18=$(echo "$count_ref_pop18 + $count_alt_pop18" | bc)
		elif [ "$majorall_pop18" == "$minorall_pop18" ]; then
			count_ref_pop18=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_alt_pop18=$(echo "NA") # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
			count_ref_alt_pop18=$(echo "NA") 
			echo "cas 2: majorall_pop18 == minorall_pop18 pour $majorall_pop18 == $minorall_pop18 à la position $position" >> $pwd_biallelic_sites_controlbugs
		elif [ "$minorall_pop18" == "N" ]; then
			echo "cas 3: minorall_pop1 == N où $majorall_pop18 != $ref à la position $position" >> $pwd_biallelic_sites_controlbugs
			if [ "$ref" == "$majorall_pop18" ]; then
				count_ref_pop18=$(awk '{print $81}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop18=$(awk '{print $117}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop18=$(echo "$count_ref_pop18 + $count_alt_pop18" | bc)
			else
				count_ref_pop18=$(awk '{print $117}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_alt_pop18=$(awk '{print $81}' $tmp00) # éliminer les cas "N" et "N" pour major et minor (faible couverture?)
				count_ref_alt_pop18=$(echo "$count_ref_pop18 + $count_alt_pop18" | bc)
			fi
		else
			count_ref_pop18=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_alt_pop18=$(echo "NA") # cas où l'allèle de référence ne serait ni le major ni le minor, devrait être très très très rare
			count_ref_alt_pop18=$(echo "NA")
			echo "autre cas à analyser: majorall_pop18=$majorall_pop18 minorall_pop18=$minorall_pop18 ref=$ref position=$position" >> $pwd_biallelic_sites_controlbugs
		fi
	else # 3eme allele retenu à ce locus alors que le site est attendu biallelic, ça veut probablement suggérer un locus monomorphe pour le major et une erreur de séquençage retenue. Néanmoins dans le doute je ne conserve pas ces sites
		count_ref_pop18=$(echo "NA") 
		count_alt_pop18=$(echo "NA")
		count_ref_alt_pop18=$(echo "NA")
		echo "cas 1: minorall_pop18 == $minorall_pop18 != $expectedmajor ou $expectedminor à la position $position" >> $pwd_biallelic_sites_controlbugs
	fi
	echo "position=$position	count_ref_pop1=$count_ref_pop1	count_alt_pop1=$count_alt_pop1	count_total_pop1=$count_ref_alt_pop1"
	if [ "$count_alt_pop1" == "NA" ] || [ "$count_alt_pop2" == "NA" ] || [ "$count_alt_pop3" == "NA" ] || [ "$count_alt_pop4" == "NA" ] || [ "$count_alt_pop5" == "NA" ] || [ "$count_alt_pop6" == "NA" ] || [ "$count_alt_pop7" == "NA" ] || [ "$count_alt_pop8" == "NA" ] || [ "$count_alt_pop9" == "NA" ] || [ "$count_alt_pop10" == "NA" ] || [ "$count_alt_pop11" == "NA" ] || [ "$count_alt_pop12" == "NA" ] || [ "$count_alt_pop13" == "NA" ] || [ "$count_alt_pop14" == "NA" ] || [ "$count_alt_pop15" == "NA" ] || [ "$count_alt_pop16" == "NA" ] || [ "$count_alt_pop17" == "NA" ] || [ "$count_alt_pop18" == "NA" ]; then
		echo "$position: no test because this site contains NA information at least for one population" >> $pwd_biallelic_sites_counts_infotest
		echo "$position: this site contains NA information at least for one population - sites were excluded for this version of the script alt_pop1= $count_alt_pop1 alt_pop2= $count_alt_pop2  alt_pop3= $count_alt_pop3 alt_pop4= $count_alt_pop4 alt_pop5= $count_alt_pop5 alt_pop6= $count_alt_pop6 alt_pop7= $count_alt_pop7 alt_pop8= $count_alt_pop8 alt_pop9= $count_alt_pop9 alt_pop10= $count_alt_pop10 alt_pop11= $count_alt_pop11 alt_pop12= $count_alt_pop12 alt_pop13= $count_alt_pop13 alt_pop14= $count_alt_pop14 alt_pop15= $count_alt_pop15 alt_pop16= $count_alt_pop16 alt_pop17= $count_alt_pop17 alt_pop18= $count_alt_pop18" >> $pwd_biallelic_sites_counts_infosites
		echo "$position	$count_ref_pop1	$count_alt_pop1	$count_ref_alt_pop1	$count_ref_pop2	$count_alt_pop2	$count_ref_alt_pop2	$count_ref_pop3	$count_alt_pop3	$count_ref_alt_pop3	$count_ref_pop4	$count_alt_pop4	$count_ref_alt_pop4	$count_ref_pop5	$count_alt_pop5	$count_ref_alt_pop5	$count_ref_pop6	$count_alt_pop6	$count_ref_alt_pop6	$count_ref_pop7	$count_alt_pop7	$count_ref_alt_pop7	$count_ref_pop8	$count_alt_pop8	$count_ref_alt_pop8	$count_ref_pop9	$count_alt_pop9	$count_ref_alt_pop9	$count_ref_pop10	$count_alt_pop10	$count_ref_alt_pop10	$count_ref_pop11	$count_alt_pop11	$count_ref_alt_pop11	$count_ref_pop12	$count_alt_pop12	$count_ref_alt_pop12	$count_ref_pop13	$count_alt_pop13	$count_ref_alt_pop13	$count_ref_pop14	$count_alt_pop14	$count_ref_alt_pop14	$count_ref_pop15	$count_alt_pop15	$count_ref_alt_pop15	$count_ref_pop16	$count_alt_pop16	$count_ref_alt_pop16	$count_ref_pop17	$count_alt_pop17	$count_ref_alt_pop17	$count_ref_pop18	$count_alt_pop18	$count_ref_alt_pop18" >> $pwd_biallelic_sites_counts_excludedNA
	else
		total_count_alt=$(echo "$count_alt_pop1" + "$count_alt_pop2" + "$count_alt_pop3" + "$count_alt_pop4" + "$count_alt_pop5" + "$count_alt_pop6" + "$count_alt_pop7" + "$count_alt_pop8" + "$count_alt_pop9" + "$count_alt_pop10" + "$count_alt_pop11" + "$count_alt_pop12" + "$count_alt_pop13" + "$count_alt_pop14" + "$count_alt_pop15" + "$count_alt_pop16" + "$count_alt_pop17" + "$count_alt_pop18" | bc)
		total_all_count=$(echo "$count_ref_alt_pop1" + "$count_ref_alt_pop2" + "$count_ref_alt_pop3" + "$count_ref_alt_pop4" + "$count_ref_alt_pop5" + "$count_ref_alt_pop6" + "$count_ref_alt_pop7" + "$count_ref_alt_pop8" + "$count_ref_alt_pop9" + "$count_ref_alt_pop10" + "$count_ref_alt_pop11" + "$count_ref_alt_pop12" + "$count_ref_alt_pop13" + "$count_ref_alt_pop14" + "$count_ref_alt_pop15" + "$count_ref_alt_pop16" + "$count_ref_alt_pop17" + "$count_ref_alt_pop18" | bc)
		cutoff=$(echo "$cutoffMAF * $total_all_count" | bc)
		cutoff2=$(echo "$total_all_count - ($cutoffMAF * $total_all_count)" | bc) # par sécurité, normalement ça ne doit pas se produire car je l'intègre déjà dans les loops
		test=$(echo "$total_count_alt > $cutoff && $total_count_alt < $cutoff2 && $cutoff > 5" | bc)
		if [ $test = '1' ]; then
			echo  "$position: test ACCEPTED cutoff=$cutoff total_count_alt=$total_count_alt total_all_count=$total_all_count" >> $pwd_biallelic_sites_counts_infotest
			echo "$position accepted car $total_count_alt (avec $total_all_count * $cutoffMAF) est supérieur à $cutoff" >> $pwd_biallelic_sites_counts_infosites	
			echo "$position	$count_ref_pop1	$count_alt_pop1	$count_ref_alt_pop1	$count_ref_pop2	$count_alt_pop2	$count_ref_alt_pop2	$count_ref_pop3	$count_alt_pop3	$count_ref_alt_pop3	$count_ref_pop4	$count_alt_pop4	$count_ref_alt_pop4	$count_ref_pop5	$count_alt_pop5	$count_ref_alt_pop5	$count_ref_pop6	$count_alt_pop6	$count_ref_alt_pop6	$count_ref_pop7	$count_alt_pop7	$count_ref_alt_pop7	$count_ref_pop8	$count_alt_pop8	$count_ref_alt_pop8	$count_ref_pop9	$count_alt_pop9	$count_ref_alt_pop9	$count_ref_pop10	$count_alt_pop10	$count_ref_alt_pop10	$count_ref_pop11	$count_alt_pop11	$count_ref_alt_pop11	$count_ref_pop12	$count_alt_pop12	$count_ref_alt_pop12	$count_ref_pop13	$count_alt_pop13	$count_ref_alt_pop13	$count_ref_pop14	$count_alt_pop14	$count_ref_alt_pop14	$count_ref_pop15	$count_alt_pop15	$count_ref_alt_pop15	$count_ref_pop16	$count_alt_pop16	$count_ref_alt_pop16	$count_ref_pop17	$count_alt_pop17	$count_ref_alt_pop17	$count_ref_pop18	$count_alt_pop18	$count_ref_alt_pop18" >> $pwd_biallelic_sites_counts
			echo "$position	$count_ref_pop1	$count_ref_alt_pop1	$count_ref_pop2	$count_ref_alt_pop2	$count_ref_pop3	$count_ref_alt_pop3	$count_ref_pop4	$count_ref_alt_pop4	$count_ref_pop5	$count_ref_alt_pop5	$count_ref_pop6	$count_ref_alt_pop6	$count_ref_pop7	$count_ref_alt_pop7	$count_ref_pop8	$count_ref_alt_pop8	$count_ref_pop9	$count_ref_alt_pop9	$count_ref_pop10	$count_ref_alt_pop10	$count_ref_pop11	$count_ref_alt_pop11	$count_ref_pop12	$count_ref_alt_pop12	$count_ref_pop13	$count_ref_alt_pop13	$count_ref_pop14	$count_ref_alt_pop14	$count_ref_pop15	$count_ref_alt_pop15	$count_ref_pop16	$count_ref_alt_pop16	$count_ref_pop17	$count_ref_alt_pop17	$count_ref_pop18	$count_ref_alt_pop18" >> $pwd_biallelic_sites_counts2n
			echo "$position	$count_alt_pop1	$count_ref_alt_pop1	$count_alt_pop2	$count_ref_alt_pop2	$count_alt_pop3	$count_ref_alt_pop3	$count_alt_pop4	$count_ref_alt_pop4	$count_alt_pop5	$count_ref_alt_pop5	$count_alt_pop6	$count_ref_alt_pop6	$count_alt_pop7	$count_ref_alt_pop7	$count_alt_pop8	$count_ref_alt_pop8	$count_alt_pop9	$count_ref_alt_pop9	$count_alt_pop10	$count_ref_alt_pop10	$count_alt_pop11	$count_ref_alt_pop11	$count_alt_pop12	$count_ref_alt_pop12	$count_alt_pop13	$count_ref_alt_pop13	$count_alt_pop14	$count_ref_alt_pop14	$count_alt_pop15	$count_ref_alt_pop15	$count_alt_pop16	$count_ref_alt_pop16	$count_alt_pop17	$count_ref_alt_pop17	$count_alt_pop18	$count_ref_alt_pop18" >> $pwd_biallelic_sites_counts2n
		elif  [ $test != '1' ]; then
			echo "$position: test REJECTED cutoff=$cutoff total_count_alt=$total_count_alt total_all_count=$total_all_count" >> $pwd_biallelic_sites_counts_infotest
			echo "$position excluded car $total_count_alt (avec $total_all_count * $cutoffMAF) est inférieur à $cutoff" >> $pwd_biallelic_sites_counts_infosites
			echo "$position	$count_ref_pop1	$count_alt_pop1	$count_ref_alt_pop1	$count_ref_pop2	$count_alt_pop2	$count_ref_alt_pop2	$count_ref_pop3	$count_alt_pop3	$count_ref_alt_pop3	$count_ref_pop4	$count_alt_pop4	$count_ref_alt_pop4	$count_ref_pop5	$count_alt_pop5	$count_ref_alt_pop5	$count_ref_pop6	$count_alt_pop6	$count_ref_alt_pop6	$count_ref_pop7	$count_alt_pop7	$count_ref_alt_pop7	$count_ref_pop8	$count_alt_pop8	$count_ref_alt_pop8	$count_ref_pop9	$count_alt_pop9	$count_ref_alt_pop9	$count_ref_pop10	$count_alt_pop10	$count_ref_alt_pop10	$count_ref_pop11	$count_alt_pop11	$count_ref_alt_pop11	$count_ref_pop12	$count_alt_pop12	$count_ref_alt_pop12	$count_ref_pop13	$count_alt_pop13	$count_ref_alt_pop13	$count_ref_pop14	$count_alt_pop14	$count_ref_alt_pop14	$count_ref_pop15	$count_alt_pop15	$count_ref_alt_pop15	$count_ref_pop16	$count_alt_pop16	$count_ref_alt_pop16	$count_ref_pop17	$count_alt_pop17	$count_ref_alt_pop17	$count_ref_pop18	$count_alt_pop18	$count_ref_alt_pop18" >> $pwd_biallelic_sites_counts_excludedMAF
		else
			echo "$position: test REJECTED (other) cutoff=$cutoff total_count_alt=$total_count_alt total_all_count=$total_all_count"  >> $pwd_biallelic_sites_counts_infotest
			echo "ELSE à $position excluded car $total_count_alt (avec $total_all_count * $cutoffMAF) pour une raison inconne cutoff = $cutoff alt_pop1= $count_alt_pop1 alt_pop2= $count_alt_pop2  alt_pop3= $count_alt_pop3 alt_pop4= $count_alt_pop4 alt_pop5= $count_alt_pop5 alt_pop6= $count_alt_pop6 alt_pop7= $count_alt_pop7 alt_pop8= $count_alt_pop8 alt_pop9= $count_alt_pop9 alt_pop10= $count_alt_pop10 alt_pop11= $count_alt_pop11 alt_pop12= $count_alt_pop12 alt_pop13= $count_alt_pop13 alt_pop14= $count_alt_pop14 alt_pop15= $count_alt_pop15 alt_pop16= $count_alt_pop16 alt_pop17= $count_alt_pop17 alt_pop18= $count_alt_pop18" >> $pwd_biallelic_sites_counts_infosites
		fi
	fi
done < $pwd_biallelic_sites_extended
