import sys
import argparse
import pysam
import random
import math
import numpy as np

#############################    ARGUMENTS     ##################################

def parse_args(args):
	ap = argparse.ArgumentParser(description = 'call consensus reads from paired-end, Duplex-UMI reads, preserving mapping information')

   # INPUT FILE
	ap.add_argument("-i", "--input_file",
			 required = True,
			 type = str,
			 help = "REQUIRED: name of input file in .bam format")

   # OUTPUT FILENAME
	ap.add_argument("-o", "--output_file",
			required = False, default = None,
			type = str,
			help = "Name of the output file with consensus reads - in .bam format \n default = \"<input_filename>_cons.bam\" ")

   # VERBOSE
	ap.add_argument("-v", "--verbose",
			action = 'store_true',
			help = "Verbose: add it as argument if you want more information on how the program is running \n Note it might take longer too \n default = False") 

   # THRESHOLD MIN MAPPING QUALITY
	ap.add_argument("-q", "--min_map_quality",
			required = False, default = 20,
			type = int,
			help = "min mapping quality of a read to use for the generation of a consensus read \n given as Phred score \n default = 20")

   # THRESHOLD MIN BASE QUALITY
	ap.add_argument("--min_base_quality",
			required = False, default = 20,
			type = int,
			help = "min base quality for a base not to be masked \n given as Phred score \n default = 20")

   # THRESHOLD MIN READS to create consensus
	ap.add_argument("--min_reads",
			required = False, default = 1,
			type = int,
			help = "min number of input reads to generate a consensus read \n default = 1")

   # THRESHOLD MAX READS to create consensus
	ap.add_argument("--max_reads",
			required = False, default = 100,
			type = int,
			help = "maximum number of input reads to generate a single-strand consensus read \n default = 100")

   # MAX BASE QUALITY ALLOWED 	
	ap.add_argument("--max_base_quality",
			required = False, default = 60,
			type = int,
			help = "all bases with a quality greater than this one will be capped to this value before calling consensus \n given as Phred score \n recommended to use in case of overestimated base qualities \n default = 60")

   # BASE QUALITY SHIFT VALUE
	ap.add_argument("--base_quality_shift",
			required = False, default = 0,
			type = int,
			help = "all input base qualities will be substracted this value before calling consensus \n given as Phred score \n recommended to use in case of overestimated base qualities \n default = 0")

   # ERROR RATE POST-LABELLING
	ap.add_argument("--error_rate_post_labeling",
			required = False, default = 0,
			type = int,
			help = "error rate for an error after the UMIs have been integrated to the fragment, but prior to sequencing (does not include the error during sequencing). Such an error could be introduced by target capture or amplification. \n given as Phred score \n recommended to use in case of having prior information about this error \n default = 0")

   # ERROR RATE PRE-LABELLING
	ap.add_argument("--error_rate_pre_labeling",
			required = False, default = 0,
			type = int,
			help = "error rate for an error before the UMIs have been integrated to the fragment. Such an error could be introduced by deamination of the base in the template DNA or oxidation of the base during library preparation.  \n given as Phred score \n recommended to use in case of having prior information about this error \n default = 0")

   # DELETION SCORE
	ap.add_argument("--deletion_score",
			required = False, default = 30,
			type = int,
			help = "quality score for a deletion \n given as a Phred score \n default = 30")

   # ABSENCE OF INSERTION SCORE
	ap.add_argument("--no_insertion_score",
			required = False, default = 30,
			type = int,
			help = "quality score for the absence of an insertion in a position where other reads have an insertion \n given as a Phred score \n default = 30")


	return (ap.parse_args(args))



################################   SMALL FUNCTIONS      ##########################################

# These functions are named "small" as they are called by other functions: the so-called "STEP FUNCTIONS"


def check_family_UMIs(family, family_code):
	'''
	Checks whether all reads from a family have the same UMIs
	UMIs are stored in "RX" tag, either as UMI1-UMI2 or UMI2-UMI1
	Raises an error if this condition is not fullfilled
	'''

	family_umi1 = family[0].get_tag("RX")
	family_umi2 = "-".join([family_umi1.split("-")[1], family_umi1.split("-")[0]])

	for read in family:
		if read.get_tag("RX") != family_umi1 and read.get_tag("RX") != family_umi2:
			print("ERROR: family", family_code, "has different UMI tags. \n Please check output file of previous step of the pipeline (fgbio GroupReadsByUmi)")
			sys.exit(1)


def check_family_rnames(family, family_code):
	'''
	Checks whether all reads from a family have the same value in rname field. e.g. same chr number
	Raises an error if this condition is not fullfilled
	'''

	# reference rname: the rname of the first read in the family
	rname = family[0].reference_id
 
	for read in family[1:]:
		if read.reference_id != rname:
			print("ERROR: family", family_code, "has difference rnames (e.g. chromosome numbers). \n Please check output file of previous step of the pipeline (fgbio GroupReadsByUmi)")
			sys.exit(1)



def split_family(family):
	'''
	Splits the input family (given as list of reads) into a list of 4 lists of reads, corresponding to each subfamily
	Returns: a list of 4 lists of reads
	The first sub-list contains reads from A1 subfamily (maps forward strand, first paired-end)
	The second sub-list contains reads from B2 subfamily (maps forward strand, second paired-end)
	The third sub-list contains reads from B1 subfamily (maps reverse strand, first paired-end)
	The fourth sub-list contains reads from A2 subfamily (maps reverse strand, second paired-end)
	'''

	family_splitted = [[], [], [], []]      # list of list of reads
						# 4 subgroups: A1, B2, B1, A2
	for read in family:
		if   not read.is_reverse  and  read.is_read1:  # read belongs to A1
			family_splitted[0].append(read)
		elif not read.is_reverse  and  read.is_read2:  # read belongs to B2
			family_splitted[1].append(read)
		elif     read.is_reverse  and  read.is_read1:  # read belongs to B1
			family_splitted[2].append(read)
		elif     read.is_reverse  and  read.is_read2:  # read belongs to A2
			family_splitted[3].append(read)

	return family_splitted


def check_number_reads(family_splitted, min_reads, max_reads, family_code, split_dict):
	'''
	Checks:
	(i) whether every subfamily in a family has enough reads to form a consensus read
	(ii) whether any subfamily has too many reads to generate a consensus in a reasonable time

	Requires:
	- family_splitted: a list of four sub-lists containing reads for each of the subfamilies
	- threshold for minimum reads
	- threshold for maximum reads
	- the code of the family, just to print a warning
	- a dictionary with the subfamilies code (which had been declared as global in main() ), just to print a warning
	- NOTE it also uses a global variable called verbose, given as argument

	Returns:
	(i) None (empty object), if any subfamily does not have enough reads
	(ii) A downsampled family_splitted for the subfamily that had too many reads. Randomly downsampled to the number max_reads given as argument
	(iii) The splitted family unmodified, if any of the above was True
	'''

	for idx, subfamily in enumerate(family_splitted):
		if len(subfamily) < min_reads:
			if verbose:
				print ("Not enough reads in family", family_code, "to form consensus reads \n")
			return None

		elif len(subfamily) > max_reads:
			if verbose:
				print ("Family", family_code, "subfamily", split_dict[idx], "has been randomly downsampled to", max_reads, "reads to form consensus read")
			family_splitted[idx] = random.sample(family_splitted[idx], max_reads)

	return family_splitted


def remove_clipping(read):
	'''
	Requires: a read object (pysam.AlignedSegment object)
	
	Removes soft and hard clips from a read
	In case of hard clips, they are removed from: CIGAR string
	In case of soft clips, they are removed from: CIGAR string, sequence, base qualities 
	There is a commented-out option to remove soft-clips from starting position, in case of 5'end soft-clips
	It assumes that soft-clips can only occur at the start/end of the sequence, but this is supposed to have been checked previously in the program, in function pass_filters()
	
	Returns: a read object (pysam.AlignedSegment object) without soft and hard clips

	Code taken from: https://github.com/ngsutils/ngsutils/blob/master/ngsutils/bam/removeclipping.py
	'''

	modified_cigar = []
	softclip_5 = 0
	softclip_3 = 0

	modified = False
	inseq = False

	### PARSE CIGAR ###
	for op, length in read.cigartuples:
		if op == 5:                     # H
			modified = True
		elif op == 4:                   # S - it assumes S can only occur at the start/end of seq
			modified = True
			if not inseq:
				softclip_5 = length
			else:
				softclip_3 = length
		else:                           # M, D or I - the remaining operators must have been filtered out before
			inseq = True
			modified_cigar.append((op, length))


	### CHANGE CIGAR ###      # in case of both S and H
	if not modified:
		return read
	read.cigartuples = modified_cigar
	original_length = len(read.query_sequence)


	### TRIM POS, SEQUENCE, SEQUENCE QUALITIES ###    # in case of S
	# COMMENT: position should be modified when trimming soft clips? for now, it's commented out

	seq  = read.query_sequence
	seqQ = read.query_qualities
	pos  = read.reference_start

	if softclip_5:
		# read.reference_start = pos + softclip_5
		pass

	if softclip_3:
		read.query_sequence  = seq  [softclip_5 : -softclip_3]
		read.query_qualities = seqQ [softclip_5 : -softclip_3]
	else:
		read.query_sequence  = seq  [softclip_5 : ]
		read.query_qualities = seqQ [softclip_5 : ]


	### ADD TRIMMING TAGS ###                       # in case of S
	newtags = []
	if softclip_5:
		newtags.append(('ZA', softclip_5))
	if softclip_3:
		newtags.append(('ZB', softclip_3))
	newtags.append(('ZC', float(softclip_5 + softclip_3) / original_length))

	read.set_tags(read.get_tags() + newtags)


	return read


def mask_low_quality_bases(read, seqQ_threshold):
	'''
	Requires: a read object (pysam.AlignedSegment object) and a Phred score value for the minimum base quality allowed
	Masks (transforms them to "N") all bases from a read having a quality smaller than the given threshold
	Qualities of masked bases are not modified (should they be transformed to 2 (minimum Phred score)?)
	Returns: a read object (pysam.AlignedSegment object) with its low-quality bases masked 	
	'''

	masked_sequence = []
	base_qualities = read.query_qualities

	for idx, quality in enumerate(read.query_qualities):
		if quality < seqQ_threshold:
			masked_sequence.append("N")
		else:
			masked_sequence.append(read.query_sequence[idx])

	read.query_sequence = "".join(masked_sequence)
	read.query_qualities = base_qualities


	return read


def trim_3prime_N(read):
	'''
	Requires: a read object (pysam.AlignedSegment object)
        Trims all the Ns located in the 3' end of the given read, modyfing the read's sequence, base qualities and CIGAR string
	Uses internal functions "expand_cigartuples()" and "compress_cigarlist()"
	Returns: a read object (pysam.AlignedSegment object)
	'''

	sequence = list(read.query_sequence)
	base_qualities = read.query_qualities

	original_seq_length = len(sequence)

	trim_3 = len(sequence)
	for nt in reversed(sequence):
		if nt == 'N' :
			trim_3 -= 1
			if verbose:
				print("3-prime N trimming has been performed on read", read.query_sequence)
		else:
			break

	sequence = sequence [: trim_3]
	base_qualities = base_qualities [: trim_3]
	
	read.query_sequence  = "".join(sequence)
	read.query_qualities = base_qualities

	original_cigar = expand_cigartuples (read.cigartuples)
	trim_cigar = len(original_cigar) - (original_seq_length - trim_3)
	modified_cigar = compress_cigarlist ( original_cigar [: trim_cigar] )
	read.cigartuples = modified_cigar 

	return read



def load_next_family (next_family):
	'''
	When a family has been read and its consensus has been generated (or it did not have enough reads to generate a consensus), 
	the read sequence from the next family, that had been already read (verb "to read"), is changed from the list next_family to the list family
	and next_family list is emptied
	also, the code of the new family is updated
	'''	
	
	family = next_family
	read = family[0]
	family_code = read.get_tag("MI").split("/")[0]
	next_family = None

	return family, family_code, next_family



def expand_cigartuples (cigartuples):
	'''
	Returns: list of CIGAR operations (corresponding to a single read), containing one operation per nucleotide in the read's sequence
	Requires: list of CIGAR tuples (corresponding to a single read)
	This is used prior to reconstructing the alignment of the reads, since the function "reconstruct_alignment" needs the CIGAR operations expanded (one number per nucleotide) not zipped in tuples
	'''
	cigarlist = []

	for operation, length in cigartuples:
		for i in range(length):
			cigarlist.append(operation)

	return cigarlist


def change_match_mismatch_operations (cigarlist):
	'''
	Returns: a list of CIGAR operations, with operations 7 (= or sequence match) and 8 (X or sequence mismatch) transformed to 1 (M or alignment match)
	The output list corresponds to a single read's CIGAR
	This is used prior to reconstructing the alignment of the reads, since the function "reconstruct_alignment" cannot handle "=" and "X" symbols.
	'''
	change = False
	
	for idx, operation in enumerate(cigarlist):
		if operation == 7 or operation == 8:
			change = True
			cigarlist[idx] = 0

	if change == True:
		print('Symbols "x" and "=" were found in CIGAR strings and were changed to "M"')

	return cigarlist


def extract_info_from_reads (list_of_reads):    
	'''
	Returns: lists of (i) positions, (ii) CIGAR tuples, (iii) sequences and (iv) base qualities in Phred scores; of all reads given as input
	This is used prior to reconstructing the alignment of the reads, to extract all the fields that the function "reconstruct_alignment" needs as arguments.
	'''

	pos_reads   = []
	cigar_reads = []
	seq_reads   = []
	qual_reads  = []
	mapQ_reads  = []

	for read in list_of_reads:
		pos_reads  .append(read.reference_start)
		cigar_reads.append(read.cigartuples)
		seq_reads  .append(read.query_sequence)
		qual_reads .append(list(read.query_qualities))
		mapQ_reads .append(read.mapping_quality)

	cigar_reads = [expand_cigartuples(x) for x in cigar_reads]
	cigar_reads = [change_match_mismatch_operations(x) for x in cigar_reads]

	seq_reads = [list(x) for x in seq_reads]
	
	
	return pos_reads, cigar_reads, seq_reads, qual_reads, mapQ_reads


def get_current_CIGAR_operations(cigar_reads, idx_cigar, number_reads):
	'''
	Returns: a list of CIGAR operations in the current position of every read
	This is done to check whether there are any insertion (op = 1) in the current position of any of the reads
	'''

	operations = []

	for x in range (number_reads):
		try:
			op = cigar_reads[x][idx_cigar[x]]
			operations.append(op)

		# an IndexError can arise in case of reads that have already reached their 3' end. 
		# if this happens, for sure there is no insertion on that position of that read
		# so an arbitrary number different than 1 is added.
		except IndexError:
			operations.append(0)

	return operations
            

def reconstruct_alignment (pos_reads, cigar_reads, seq_reads, qual_reads):
	'''
	Returns: (i) list of aligned sequences, (ii) list of aligned base qualities (Phred scores)

	Requires: lists of (i) positions, (ii) CIGAR operations, (iii) sequences and (iv) base qualities(Phred score) of each read that will be aligned.

	This function reconstructs the alignment of the reads given, based on their CIGAR strings.
	This is done to properly generate a consensus read without losing any mapping information.
	The alignment is constructed iteratively, starting on the lowest position to which any of the reads map, until the greatest position to which any of the reads map.
	It only considers CIGAR strings with (mis)matches (M), insertions (I) and deletions (D). 
	It assumes the rest possible operations (e.g. soft/hard clips) have been filtered/trimmed.
	Note that "X" and "=" operations are not considered here either and they have to be transformed to "M" before i.e. use function "change_match_mismatch_operations".
	FUTURE MODIFICATION: allow for "X" and "=" operations.

	First, it checks for insertions in the corresponding position of all reads. 
	If there are, it adds the nucleotide letter IN LOWERCASE to the aligned sequence of the read that had the insertion. 
	The reads that do not have an insertion in that position, will be added a "+" sign.
	This arbitrary method allows for a correct assesment of insertions when generating a consensus sequence and a consensus CIGAR sequence.

	If there are no insertions, it will add the nucleotide letter to the aligned sequence, in case of a match;
	or a '-' sign, in case of a deletion.

	Loop will move forward in the read's CIGAR string in case of: match, deletion and insertion ONLY in the read that had the insertion.
	Loop will move forward in the read's sequence in case of: match and insertion of ONLY in the read that had the insertion.
	Loop will not move forward in the read's CIGAR string in case of: insertion in any other read and in positions where the read does not map.
	Loop will not move forward in the read's sequence in case of: deletion, insertion in any other read and in positions where the read does not map.
	'''

	minimum_pos = min(pos_reads) # lowest position to which any of the reads map
	maximum_pos = max([pos + len(seq) for pos, seq in zip(pos_reads, seq_reads)]) # highest position to which any of the reads maps

	number_reads = len(pos_reads)

	# Indexes (init to 0) that will be updated on every iteration
	# they point to the next position of every read that should be added to the alignment
	idx_cigar = [0] * number_reads
	idx_seq   = [0] * number_reads

	# Aligned sequences and aligned qualities
	aligned_seq  = [ [] for i in range(number_reads) ]
	aligned_qual = [ [] for i in range(number_reads) ]


	for pos in range(minimum_pos, maximum_pos):

		# First check whether there are insertions in corresponding positions
		current_operations = get_current_CIGAR_operations(cigar_reads, idx_cigar, number_reads)
		# There are INSERTIONS in corresponding positions
		if 1 in current_operations:
			for read_ in range(number_reads):

				# The insertion is in this read
				if current_operations[read_] == 1:
					# add nucleotide but lowercase
					# move indexes forward
					aligned_seq [read_].append(seq_reads [read_][idx_seq[read_]] .lower())
					aligned_qual[read_].append(qual_reads[read_][idx_seq[read_]] )
					idx_cigar[read_] += 1
					idx_seq  [read_] += 1


				# The insertion was on other read  
				else:
					# add '+' to seq
					# add '+' to base qualities
					# do not move any index forward
					aligned_seq [read_].append('+')
					aligned_qual[read_].append('+')

			continue
    
		# There are no insertions in corresponding positions
		for read_ in range(number_reads): 

			# Read does not map to the current position
			# 5' end has not been reached yet
			if pos < pos_reads[read_]:
				# add N to seq
				# add 2 (minimum Phred score) to base qualities
				aligned_seq [read_].append('N')
				aligned_qual[read_].append(2)


			# Read maps to the current position
			elif idx_seq[read_] < len(seq_reads[read_]):

				# DELETION on current read and current position
				if cigar_reads[read_][idx_cigar[read_]] == 2:
					# add '-' to seq
					# add '-' to base qualities
					# move only cigar index forward
					# do not move seq index                    
					aligned_seq [read_].append('-')
					aligned_qual[read_].append('-')
					idx_cigar [read_] += 1


				# MATCH on current read and current position
				else:
					# add nucleotide to seq
					# add score to base qualities
					# move indexes forward
					aligned_seq [read_].append(seq_reads [read_][idx_seq[read_]])
					aligned_qual[read_].append(qual_reads[read_][idx_seq[read_]])
					idx_cigar[read_] += 1
					idx_seq  [read_] += 1


			# Read does not map to this position
			# 3' end has been already reached
			else: 
				# add N to seq
				# add 2 (min Phred score) to base qualities
				aligned_seq [read_].append('N')
				aligned_qual[read_].append(2)


	return aligned_seq, aligned_qual, minimum_pos


def most_likely_nucleotide(nucleotides, prob_error, threshold_posterior_prob):
	'''
	Calculates likelihoods of each nucleotide (A,T,C,G) as well as of a deletion and an insertion given a list of nucleotides corresponding to the nucleotides of every read in a given position.
	Based on these likelihoods, it returns the most likely nucleotide/INDEL.
	It is used to generate the consensus sequence, to determine the most likely nucleotide/INDEL on the given position

	Requires:
	(i) list of nucleotides corresponding to the bases from a single specific position of every input read used to form the consensus
	possible values are: 
		- 'A', 'T', 'C', 'G', 'N',
		- 'a', 't', 'c', 'g', 'n' (lowercase corresponds to insertions)
		- '+' (lack of insertion in that read, but insertion in other read in that specific position) 
		- '-' (deletion)
	(ii) base qualities for each input nucleotide (i) transformed to probabilities of an error during sequencing (this probability might also include the probability of an error before sequencing, see call_consensus function)
	(iii) a threshold for the posterior probability of the here-calculated most likely nucleotide/INDEL. Below this threshold, the most likely nucleotide will be masked.

	Returns:
	(i) single nucleotide/insertion/deletion, corresponding to the nucleotide/indel with the highest posterior probability among all possible values. 
	Possible output values are the same as the possible input values for argument (i).
	(ii) posterior probability of the most likely nucleotide

	Likelihoods and posteriors are calculated based on fulcrumgenomics benchmark, but also including INDELs. For more information on the mathematics behind this implementation, please see documentation or visit https://github.com/fulcrumgenomics/fgbio/wiki/Calling-Consensus-Reads.

	Note that the likelihood of "a" ("A" on insertion) is considered together with the likelihood of "A" ("A" on sequence match). Whether the output most likely nucleotide is lowercase (insertion) or uppercase (sequence match) will depend on whether the given list of nucleotides only contained insertions (lowercase nucletodides, and '+' in the reads where the insertion was not present) or only contained deletions/sequence matches.
	Note that likelihood of "N" is not considered here, since "N" is not a biological base or operation that can be observed in a sequence.
	However, "N" is one of the possible outputs of the function, only if the posterior probability of the most likely nucleotide/INDEL is below the given threshold.
	'''

	###   0. Check input nucleotide values are correct   ###

	possible_values = ['A', 'T', 'C', 'G', 'N', 'a', 't', 'c', 'g', 'n', '+', '-']

	if bool ( [nt for nt in nucleotides if (nt not in possible_values) ]):
		print ("nucleotides = ", nucleotides)
		print ("\n ERROR: input values to nucleotide/INDEL list in function call_consensus() and subfunction most_likely_nucleotide() are different than the allowed ones. \n Please check documentation") 
		sys.exit(1)


	###   1. Compute likelihoods of each nucleotide/INDEL   ###

	likelihoods = [1, 1, 1, 1, 1, 1] 				 # One position for each nucleotide / INDEL
	nt_dict = {'A':0, 'T':1, 'C':2, 'G':3, '+':4, '-':5, 'N':6}	 # Index corresponding to each nucleotide's likelihood. Note that the dictionary has 7 entries and the list "likelihoods" has only length 6, since the last entry in the dictionary (N) is not considered in the likelihoods
	nt_dict_inverse = {0:'A', 1:'T', 2:'C', 3:'G', 4:'+', 5:'-', 6:'N'}

	for nt, prob_err in zip(nucleotides, prob_error):
		nt_idx = nt_dict [nt.upper()]
		for i in range(len(likelihoods)):
			if i == nt_idx:
				likelihoods[i] *= (1 - prob_err)
			else:
				likelihoods[i] *= prob_err/(len(likelihoods) - 1)


	###   2. Compute posterior probabilities of each nucleotide/INDEL   ###

	# this is the probability of observing a nucleotide, not the probability of observing a sequencing error
	posteriors = []
	for i in range (len(likelihoods)):
		posteriors.append( likelihoods[i] / np.sum(likelihoods) ) 
	

	###   3. Determine most likely nucleotide/INDEL   ###

	most_likely_nt = nt_dict_inverse[np.argmax(posteriors)].lower() if '+' in nucleotides else nt_dict_inverse[np.argmax(posteriors)]
	posterior_prob = np.max(posteriors)
	
	# if posterior probability of most likely NT < threshold --> mask the nucleotide (return an "N")
	if posterior_prob < threshold_posterior_prob:
		most_likely_nt = "N".lower() if '+' in nucleotides else "N"


	return most_likely_nt, posterior_prob



def call_consensus(seq, qual, base_quality_shift, max_base_quality, error_rate_post_labeling, error_rate_pre_labeling, seqQ_threshold, deletion_score, no_insertion_score):
	'''
	Generates consensus sequence and consensus base qualities, given multiple reads' sequences and base qualities.

	Requires: 
	(i) List of sequences of all input reads, aligned according to their CIGAR strings and starting positions (see recontruct_alignemnt function() for more details)
	(ii) List of base qualites of all input reads, aligned according to their CIGAR strings and starting positions (see recontruct_alignemnt function() for more details)
	(iii) Base quality shift: all input base qualities will be substracted this value (default = 0) before generating consensus. This option is included (default = 0, so by default this does not change anything) in case of overestimated base qualities.
	(iv) Max base quality allowed: all input base qualities will be capped to this value (default = 60) before generating consensus. This option is included (default = 60, so by default this does not change anything) in case of overestimated base qualities. 
	(v) Error rate post labeling: error rate for an error after the UMIs have been integrated to the fragment, but prior to sequencing (does not include the error during sequencing). Such an error could be introduced by target capture or amplification. This option is included (defalt = 0, so by default this does not change anything) in case of having prior information about this error
	(vi) Error rate pre labeling: error rate for an error before the UMIs have been integrated to the fragment. Such an error could be introduced by deamination of the base in the template DNA or oxidation of the base during library preparation. This option is included (default = 0, so by default this not change anything) in case of having prior information about this error
	(vii) Base quality threshold to mask a base when its quality is lower than this value.
	(viii) Base quality score in case of a deletion. This quality score would the same for every deletion found in the sequences.
	(ix) Base quality score in case of an "absence of insertion" (where there are other reads that have an insertion in that position, those reads which do not have an insertion are considered to have an "absence of insertion" in that position, represented by a '+'). This quality score would be the same for every "absence of insertion" found in the sequences.
	
	Returns:
	(i) Consensus sequence in a list format. It is aligned following the same principles as the input sequences (see reconstruct_alignment() function)
	(ii) Consensus base qualities in a list format. They are also aligned.

	Consensus sequence and base qualities are calculated based on fulcrumgenomics benchmark, but also including INDELs. For more information on the mathematics behind this implementation, please see documentation or visit https://github.com/fulcrumgenomics/fgbio/wiki/Calling-Duplex-Consensus-Reads

	Firstly, base qualities are adjusted and transformed to probabilities of sequencing error (including errors after integrating the UMIs)
	Secondly, position by position, the most likely nucleotide or INDEL (deletion, insertion or absence of insertion (when in the current position there are insertions in other reads but not in the current read)) is calculated (function most_likely_nucleotide() is called), as well as the posterior probability of having that nucleotide/INDEL. 
	These nucleotide/INDELs will be added to the consensus sequence list
	Thirdly, consensus posterior probabilities are adjusted and back-transformed to Phred scores (errors prior to integrating the UMIs are included at this step). 
	These Phred scores will be added to the consensus base qualities list
	'''

	# Create new lists

	prob_error_post_labeling = []
	posterior_prob = []
	for i in range (len(qual)):
		prob_error_post_labeling.append (qual[i][:])
	cons_seq = []
	cons_qual = []


	####    1. ADJUST INPUT BASE QUALITIES    ####

	for i in range (len(qual)):
		for j in range (len(qual[i])):
			if qual [i][j] == '+':
				prob_seq_error = 10 ** ( -no_insertion_score / 10 )
				prob_error_post_labeling [i][j] = error_rate_post_labeling*(1-prob_seq_error) + (1-error_rate_post_labeling)*prob_seq_error + error_rate_post_labeling*prob_seq_error*4/5
			elif qual [i][j] == '-':
				prob_seq_error = 10 ** ( -deletion_score / 10 )
				prob_error_post_labeling [i][j] = error_rate_post_labeling*(1-prob_seq_error) + (1-error_rate_post_labeling)*prob_seq_error + error_rate_post_labeling*prob_seq_error*4/5
			else:
				qual_adjusted = min (qual [i][j] - base_quality_shift, max_base_quality)
				prob_seq_error = 10 ** ( -qual_adjusted / 10 )
				prob_error_post_labeling [i][j] = error_rate_post_labeling*(1-prob_seq_error) + (1-error_rate_post_labeling)*prob_seq_error + error_rate_post_labeling*prob_seq_error*4/5

	# Threshold probability to mask a nucleotide  
	threshold_prob_error = 10 ** ( -seqQ_threshold / 10)
	threshold_posterior_prob = 1 - threshold_prob_error


	####    2. COMPUTE LIKELIHOODS OF NUCLEOTIDES AND INDELS    ####

	seq  = np.array(seq)
	prob_error_post_labeling = np.array(prob_error_post_labeling).astype(np.float)      # equivalent to qualities

	for nucleotides_position, prob_error_position in zip(np.transpose(seq), np.transpose(prob_error_post_labeling)):

		nt, prob = most_likely_nucleotide(nucleotides_position, 
						  prob_error_position,
						  threshold_posterior_prob)
		cons_seq.append(nt)
		posterior_prob.append(prob)


	####    3. ADJUST CONSENSUS BASE QUALITIES    ####

	for i in range (len(posterior_prob)):
		cons_prob_error_post_labeling = 1 - posterior_prob[i]
		cons_prob_error_pre_and_post_labeling = error_rate_pre_labeling*(1-cons_prob_error_post_labeling) + (1-error_rate_post_labeling)*cons_prob_error_post_labeling + error_rate_pre_labeling*cons_prob_error_post_labeling*4/5
		
		try:
			q = int( round (-10 * math.log10 (cons_prob_error_pre_and_post_labeling), 0))
			cons_qual.append(min(q, max_base_quality))
		
		# prob_error is 0 - logarithm cannot be calculated
		except ValueError:
			cons_qual.append(max_base_quality)


	return cons_seq, cons_qual



def compress_cigarlist (cigarlist):
	'''
	Opposite function of expand_cigartuples()
	Requires: list of CIGAR operations (corresponding to a single read), containing one operation per nucleotide in the read's sequence
        Returns: list of CIGAR tuples (corresponding to a single read)
	This is used after generating the consensus CIGAR, whose initial format is a list of CIGAR operations.
	The consensus CIGAR's desired format, which is returned by this function, is a list of CIGAR tuples
	'''
	cigartuples = []
	
	previous_op = None
	
	for op in cigarlist:
		if previous_op == None:
			cigartuples.append([op, 1])
			previous_op = op
		elif op != previous_op:
			cigartuples[-1] = tuple(cigartuples[-1])
			cigartuples.append([op, 1])
			previous_op = op
		else:
			cigartuples[-1][1] += 1
			previous_op = op
	
	cigartuples[-1] = tuple(cigartuples[-1])
	
	return cigartuples


def adjust_consensus_fields (cons_seq, cons_qual, minimum_pos):
	'''
	Adjusts some fields from a consensus read that have been calculated recently.
	The purpose is to obtain the desired format of these fields to include them in a read (pysam.AlignedSegment) object
	It is called by the step function "make_consensus_read()"
	
	Requires: 
	(i) consensus sequence as a list of nucleotides, including '+' (absence of insertion), '-' (deletion) and 'N' (masked bases)
	(ii) consensus base qualities as a list of Phred scores, including '+' (absence of insertion), '-' (deletion) and 'N' (masked bases) 
	(iii) lowest position that any of the reads used to generate the consensus maps to
	(iv) FUTURE: list of mapping qualities of all the reads used to generate the consensus may be required in the future
		now it is not required as the calculation of the consensus mapping quality is done in another function
	
	Returns:
	(i) adjusted consensus sequence, as a string, including only 'ATCGN'
	(ii) adjusted consensus base qualities, as a list, including only Phred scores
	(iii) consensus CIGAR tuples, as a list of tuples containing CIGAR operations and the length of these operations
	(iv) starting position of the consensus read, as an integer
	(v) FUTURE: consensus mapping quality may be returned here in the future
		if the code to calculate consensus mapping quality is included here
		now it is done separately in another function
	'''

	###   1. Trim sequence and qualities where there are N in the ends   ###
	# 5' end
	trim_5 = 0
	for nt in cons_seq:
		if nt == 'N' :
			trim_5 += 1
		else:
			break
	# 3' end
	trim_3 = len(cons_seq)
	for nt in reversed(cons_seq):
		if nt == 'N':
			trim_3 -= 1
		else:
			break
	cons_seq  = cons_seq [trim_5 : trim_3]
	cons_qual = cons_qual[trim_5 : trim_3]



	###    2. Adjust consensus starting position   ###
	# starting_position = minimum_position + number of 5' N
	cons_pos = minimum_pos + trim_5


	###    3. Adjust consensus CIGAR string   ###
	cons_cigar = []
	next_nt_checked = False

	for idx, nt in enumerate(cons_seq):

		# Only occurs in case of an insertion and a deletion continuously
		if next_nt_checked == True:
			next_nt_checked = False
			continue
		
		# Insertion
		elif nt.islower():
			# check for deletion after insertion
			# since deletion + insertion = alignment match
			try:
				next_nt = cons_seq[idx + 1]
				# there is a deletion after the insertion, therefore alignment match
				if next_nt == '-':
					cons_cigar.append(0)
					next_nt_checked = True
				# just a normal insertion
				else:
					cons_cigar.append(1)
			except IndexError:
				cons_cigar.append(1)


		# Absence of insertion ('+' sign)
		elif nt == '+':
			pass


		# Deletion
		elif nt == '-':
			# check for insertion after deletion
			# since deletion + insertion = alignment match
			try:
				next_nt = cons_seq[idx + 1]
				# there is an insertion after the deletion, therefore alignment match
				if next_nt.islower():
					cons_cigar.append(0)
					next_nt_checked = True
				# just a normal deletion
				else:
					cons_cigar.append(2)
			except IndexError:
				cons_cigar.append(2)

		# Alignment match
		else:
			cons_cigar.append(0)


	# Turn the list of CIGAR operations into a list of CIGAR tuples 
	cons_cigar = compress_cigarlist(cons_cigar)



	###    4. Adjust sequence and base qualities   ###
	
	cons_seq  = np.array(cons_seq)
	cons_qual = np.array(cons_qual)

	# Remove '+', '-' from both seq and qual
	idx_delete = np.where((cons_seq == '+') | (cons_seq == '-'))
	cons_seq  = np.delete(cons_seq,  idx_delete)
	cons_qual = np.delete(cons_qual, idx_delete)

	# Only for sequence: change small to capital letters 
	#		     and transform to string
	cons_seq = [nt.upper() for nt in cons_seq]
	cons_seq = ''.join(cons_seq)

	# Transform base qualities back to a list
	cons_qual = list(cons_qual)

		
	return cons_seq, cons_qual, cons_cigar, cons_pos


def calculate_consensus_mapping_quality (mapQ_reads):
	'''
	Calculates the mapping quality of the consensus read based on the mapping qualities of the reads used to generate the consensus
	Requires: list of mapping qualities of the reads used to generate the consensus read
	Returns: the consensus read's mapping quality

	THIS FUNCTION NEEDS IMPROVEMENT.
	Now the consensus mapping quality is calculated based on the average mapping quality of all the reads, 
	but a more precise benchmark should be thought and implemented to get more realistic and meaningful results.
	'''

	mapQ_reads = np.array (mapQ_reads)

	cons_mapQ = np.mean (mapQ_reads)

	return cons_mapQ


def get_consensus_id (read, method):
	'''
	Returns a string which will be used as the ID/query_name of the consensus read
	
	This function can be modified and customized in any desired way

	At the moment, it uses one of the reads used to make the consensus as a reference to obtain the family code and subfamily/pair-end type
	And depending on the consensus type, passed in the argument named "method" (either "single_strand" or "double_strand"), the ID will have one format or another

        In case of the having a double_strand consensus ("double_strand" method), the following criterion is used:
	* consensus reads mapping the forward strand are named as "paired-end-1"
	* consensus reads mapping the reverse strand are named as "paired-end-2"
	Note this criterion is arbitrary, and might be changed and customized as desired
	'''

	family_code = read.get_tag("MI").split("/")[0]

	if method == "single_strand":
		if not read.is_reverse and read.is_read1:
			subfamily = 'A1'
		elif read.is_reverse and read.is_read2:
			subfamily = 'A2'
		elif read.is_reverse and read.is_read1:
			subfamily = 'B1'
		else:
			subfamily = 'B2'
	
		id = 'consensus_family' + str(family_code) + '_' + subfamily
		return id
		
	elif method == "double_strand":
		if not read.is_reverse:
			paired_end = '1'
		else:
			paired_end = '2'

		id = 'consensus_family' + str(family_code) + '_paired-end' + paired_end
		return id

	else:
		print ('No IDs were assigned to the consensus reads as the method was not specified. \n Please specify whether method is "single_strand" for single strand consensus calling or "double_strand" for double strand consensus calling')
		return id


def get_consensus_flag (read, method):
	'''
	Returns an integer which will be used as the flag of the consensus read
	
	This function can be modified and customized in any desired way - mostly if families are not grouped by fgbio's pipeline

	At the moment, it uses one of the reads used to make the consensus as a reference
	It is assumed that all reads used to form the consensus:
	(i) in case of single_strand method, they all have the same flag
	(ii) in case of double_strand method, they all map to the same strand (either forward or reverse)
	If this assumption is not fulfilled, this function should be changed to add a correct flag to the consensus read
	If the script is used in its original version, this assumption is fulfilled

	Depending on the consensus type, passed in the argument named "method" (either "single_strand" or "double_strand"), the flag might have different values

	In case of the having a double_strand consensus ("double_strand" method), the following criterion is used:
	* consensus reads mapping the forward strand are assigned the bit corresponding to "paired-end-1"
	* consensus reads mapping the reverse strand are assigned the bit corresponding to "paired-end-2"
	Note this criterion is arbitrary, and might be changed and customized as desired
	'''

	if method == "single_strand":
		flag = read.flag
	elif method == "double_strand":
		if not read.is_reverse:
			flag = 99
		else:
			flag = 147
	else:
		print ('No flags were assigned to the consensus reads as the method was not specified. \n Please specify whether method is "single_strand" for single strand consensus calling or "double strand" for double strand consensus calling')
		flag = None

	return flag

def calculate_depth_and_errors (aligned_seq, aligned_cons_seq):
	'''
	Calculates consensus read depth and errors to add is as tags in .bam file
	This is done based on the way fgbio adds tags to consensus reads
	
	Requires: 
	(i) aligned_seq: list of sequences from the reads used to generate the consensus file, aligned according to their CIGAR string (see reconstruct_alignment() function)
	(ii) aligned_cons_seq : sequence of the consensus read, aligned according to their CIGAR string (output of call_consensus() function)

	Returns:
	(i) d: consensus depth - count of bases from raw reads contributing to generate the consensus read at each position
	(ii) D: maximum consensus depth - the maximum depth of raw reads at any point in the consensus read
	(iii) M: minimum consensus depth - the minimum depth of raw reads at any point in the consensus read
	(iv) e: consensus errors - the count of bases from raw reads disagreeing with the final consensus base
	(v) E: consensus error rate - the fraction of bases in raw reads disagreeing with the final consensus bases

	For more information see: http://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html

	Note that:
	- Positions where there are deletions ('-') in consensus sequence are considered, but un-commenting the corresponding lines, positions where there are deletions can be excluded
	- Positions where there are absence of insertion ('+') in consensus sequence are excluded
	- If consensus base is N, any base from raw reads on that position not having an N will be considered as error (this might be modified and improved in the future)
	- Positions in raw reads having an N are also considered as errors (if consensus base is not N), but un-commenting the corresponding line, they can be excluded to calculate the errors
	- If depth = 0 in a specific position (all masked bases), error rate (for that position) = 1, increasing considerably the average error rate (E)
	''' 

	aligned_seq = np.array(aligned_seq)
	aligned_cons_seq = np.array(aligned_cons_seq)

	# Calculate consensus depth

	condition = (aligned_seq != 'N') & (aligned_seq != 'n') & (aligned_seq != '+') # & (aligned_seq != '-')
	d = np.sum (condition, axis=0)	 # consensus depth for each position
	d = np.delete(d, np.where(aligned_cons_seq == '+'))	 # delete positions having an "absence of insertion" in consensus sequence
	# d = np.delete(d, np.where(aligned_cons_seq == '-')	 # delete positions having a deletion in consensus sequence
	D = max(d)	# consensus maximum depth
	M = min(d)      # consensus minimum depth


	# Calculate consensus errors
	
	condition = (aligned_seq != aligned_cons_seq) # & (aligned_seq != "N") & (aligned_seq != "n")
	e = np.sum(condition, axis=0)	# consensus errors for each position
	e = np.delete(e, np.where(aligned_cons_seq == '+'))	# delete positions having an "absence of insertion" in consensus sequence
	# e = np.delete(e, np.where(aligned_cons_seq == '-')	# delete positions having a deletion in consensus sequence
	with np.errstate(all = 'ignore'):
		E = e/d
		E [d == 0] = 1
		E = np.round( np.mean(E), 3 )	# consensus error rate  
	

	return list(d), int(D), int(M), list(e), float(E)



def add_tags (method, list_of_reads, aligned_seq, aligned_cons_seq, mapQ_reads):
	'''
	Returns consensus read tags in pysam format (tuple of tuples)

	Requires:
	(i) method: either "single_strand" for single-strand consensus calling or "double_strand" for double-strand consensus calling
	(ii) list of reads (list of pysam.AlignedSegment objects) used to generate the consensus
	(iii) list of aligned sequences (see function reconstruct_alignment()) from the reads used to generate the consensus
	(iv) aligned consensus sequence (output of function call_consensus())
	(v) list of mapping qualities from the reads used to generate the consensus

	Returns: tuple of tuples containing the following tags (one tag per tuple):
	- If method is single_strand ("s" in tags stands for "single"): 
		- Family (MI tag)
		- UMI (RX tag)
		- List of mapping qualities from input reads (sQ tag)
		- (sd, sD, sM, se, sE tags) fgbio-specific tags referring to the depth and errors of the consensus. For more information, see calculate_depth_and_errors() function.

	- If method is double_strand:
		- Family (MI tag)
		- UMI (RX tag). Note that the order of the two UMI tags (there is one tag per sequence end) is arbitrary, since the double-strand consensus sequence is formed by single-strand consensus reads which have the same two UMIs but in different order/different ends
		
		The following tags are "repeated" two or three times: "a" refers to the first single-strand consensus read, "b" to the second single-strand consensus read, and "c" to the double-strand consensus read:
		- List of mapping qualities from the input reads (a/b/cQ tags)
		- Single-strand consensus sequences (ac, bc tags)
		- Single-strand consensus base qualities (aq, bq tags)
		- (a/b/cd, a/b/cD, a/b/cM, a/b/ce, a/b/cE tags) fgbio-specific tags referring to the depth and errors of the consensus. For more information, see calculate_depth_and_errors() function.	 
	'''

	if method == "single_strand":

		# Family
		mi = ("MI", list_of_reads[0].get_tag("MI"))

		# UMI
		rx = ("RX", list_of_reads[0].get_tag("RX"))

		# List of mapping qualities
		sQ = ("sQ", list(mapQ_reads))

		# Depth and errors
		d, D, M, e, E = calculate_depth_and_errors (aligned_seq, aligned_cons_seq)
		sd = ("sd", str(d))
		sD = ("sD", D)
		sM = ("sM", M)
		se = ("se", str(e))
		sE = ("sE", E)

		tags = (mi, rx, sQ, sd, sD, sM, se, sE)


	elif method == "double_strand":

		# Family - only number
		mi = ("MI", list_of_reads[0].get_tag("MI").split("/")[0] )

		# UMI (order of UMI1-UMI2 is arbitrary)
		rx = ("RX", list_of_reads[0].get_tag("RX"))

		# List of mapping qualities
		aQ = ("aQ", list(list_of_reads[0].get_tag("sQ")))
		bQ = ("bQ", list(list_of_reads[1].get_tag("sQ")))
		cQ = ("cQ", list(mapQ_reads))

		# Single-strand consensus sequences 
		ac = ("ac", list_of_reads[0].query_sequence)
		bc = ("bc", list_of_reads[1].query_sequence)

		# Single-strand consensus qualities - ASCII encoded
		aq = ("aq", ''.join(map(lambda x: chr( x+33 ), list_of_reads[0].query_qualities)))
		bq = ("bq", ''.join(map(lambda x: chr( x+33 ), list_of_reads[1].query_qualities)))

		# Depth and errors
		d, D, M, e, E = calculate_depth_and_errors (aligned_seq, aligned_cons_seq)
		# Consensus depth
		ad = ("ad", list_of_reads[0].get_tag("sd"))
		bd = ("bd", list_of_reads[1].get_tag("sd"))
		cd = ("cd", str(d))
		# Max consensus depth
		aD = ("aD", list_of_reads[0].get_tag("sD"))
		bD = ("bD", list_of_reads[1].get_tag("sD"))
		cD = ("cD", D)
		# Min consensus depth
		aM = ("aM", list_of_reads[0].get_tag("sM"))
		bM = ("bM", list_of_reads[1].get_tag("sM"))
		cM = ("cM", M)
		# Consensus errors
		ae = ("ae", list_of_reads[0].get_tag("se"))
		be = ("be", list_of_reads[1].get_tag("se"))
		ce = ("ce", str(e))
		# Consensus error rate
		aE = ("aE", list_of_reads[0].get_tag("sE"))
		bE = ("bE", list_of_reads[1].get_tag("sE"))
		cE = ("cE", E)

		tags = (mi, rx, aQ, bQ, cQ, ad, bd, cd, aD, bD, cD, aM, bM, cM, ae, be, ce, aE, bE, cE, ac, bc, aq, bq)

	else:
		print ('No tags we added to the consensus reads as the method was not specified. \n Please specify whether method is "single_strand" for single strand consensus calling or "double_strand" for double strand consensus calling')
		tags = None

	return tags


#################################################    STEP FUNCTIONS    ########################################################

# Each of these functions represent one step in these pipeline, and they are all called in main()
# These functions call, at the same time, the so-called "SMALL FUNCTIONS" seen above


def pass_filters(read):
	'''
	STEP 1 OF PIPELINE
	Checks whether a read passes certain quality filters (e.g. being mapped, properly paired, passing QC filters, having a mapping quality greater than a threshold, being primary alignment...) 
	controls for possible errors in the read format (e.g. absence of MI and RX tags, presence of P/N/B/* symbols in CIGAR strings, and presence of softclipped bases in the middle of the sequence)
	Requires: read object (pysam.AlignedSegment object)
	Returns: 
	* True, if read passes all filters
	* False, if read does not pass all filters
	* Raises an error and exits the execution if there is any error in the expected format of the read	
	'''

	# Control possible errors in read format
	try:
		read.get_tag("MI")
	except:
		print('ERROR: family code tag (MI) not found in file')
		sys.exit(1)
	try:
		read.get_tag("RX")
	except:
		print('ERROR: family code tag (RX) not found in file')
		sys.exit(1)

	if bool ([ x for x in ['P', 'N', 'B', '*'] if (x in read.cigarstring) ]):
		print("ERROR: unexpected symbols (P, N, B, *) were found in CIGAR strings.")
		sys.exit(1)

	if bool ([ x for x in read.cigartuples[1:-1] if (x[0] == 4) ]):
		print("ERROR: softclips (S) found in the middle of the read.")
		sys.exit(1)


	# Pass read filters

	filters = [read.is_paired,
		   read.is_proper_pair,
		   not read.is_unmapped, 
		   not read.mate_is_unmapped,
		   not read.is_supplementary,
		   not read.is_qcfail,
		   read.mapping_quality >= mapQ_threshold]

	if all(filters):
		return True
	else:
		return False



def add_read_to_family(read, family, family_code):
	'''
	STEP 2 OF PIPELINE
	Checks whether the read belongs to the currently-loaded family.
	If so, it adds it to the current family (list of reads)
	Otherwise, it adds it to a new family, named next_family

	Requires: 
	(i) read object (pysam.AlignedSegment object)
	(ii) list of reads belonging to the currently-loaded family. It can be empty (None) if any family has been loaded yet (case of the first read in file)
	(iii) integer representing the code of the currently-loaded family. It can be empty (None) if any family has been loaded yet (case of the first read in file) 

	Returns:
	(i) the given family (list of reads), with the new read added in case of belonging to that family. If it belongs to a new family, the exact same list it was given is returned
	(ii) the family code, exactly as it was given (unless it was given as None, then it will return the code of the family the given read belongs to)
	(iii) a list named next_family, either: empty, when the input read belongs to currently-loaded family; or containing the input read, when this read belongs to a different, new family
	'''

	if family == None:                  			# only 1st read in file
		family = [read]
		family_code = read.get_tag("MI").split("/")[0]
		next_family = None

	
	elif read.get_tag("MI").split("/")[0] == family_code:   # read belongs to same family
		family.append(read)
		next_family = None
		
	else:							# read belongs to next family, current family is completed
		next_family = [read]


	return family, family_code, next_family




def preprocess_family(family, family_code):
	'''
	STEP 3 OF PIPELINE
	Once a family (list of reads) has been completely loaded in memory, it will be processed to form a consensus
	A family has been completely loaded when the variable next_family is not empty.
	When this occurs, this function is called to proccess the currently-loaded family

	Requires:
	(i) family: a list of reads belonging to the same family
	(ii) an integer representing the code number of the input family

	Sub-steps:
	1) Control for possible errors in the format of the family (e.g. different UMI tags, different rnames). If there are mistakes, it raises an error and exits the execution.
	2) Splits the family in 4 groups, depending on whether the read maps to the forward or reverse strand and whether it is first or second paired-end
	3) Checks whether the number of reads in all 4 groups is greater than a threshold to form a consensus. If this is the case, it exits the function, returning None.
	4) Removes softclipped and hardclipped bases from reads
	5) Masks bases (--> "N") of reads with a base quality lower than a given threshold (declared as global)
	6) Trims masked bases ("N") in 3' end of reads
	7) NOT IMPLEMENTED YET, BUT POSSIBLE FUTURE OPTION: trim low-quality 3' ends (right now, only trimming of 3' masked bases) 

	Returns: a list of four lists of reads (named family_splitted), containing the reads from the given family but separated into the four subgroups (mentioned above).
	If any of the four groups did not contain enough reads (see sub-step 3), then an empty object (None) is returned. 
	'''

	# Control possible errors in family

	check_family_UMIs  (family, family_code)
	check_family_rnames(family, family_code)


	# Split family

	family_splitted = split_family(family)

	
	# Check whether number of reads are enough/too many to form consensus
	# returns "None" if number of reads in ANY subfamily is not enough
	# returns a randomly downsampled list of reads if number of reads is too many

	family_splitted = check_number_reads(family_splitted, min_reads_for_consensus, max_reads_for_consensus, family_code, split_dict)

	if family_splitted == None:
		return None


	# Preprocess read by read: 
	# remove soft and hard-clipping and 
	# mask low-quality bases
	# trim masked bases ("N") on 3' end
 
	for idx, subfamily in enumerate (family_splitted):
		for subidx, read in enumerate (subfamily):

			read = remove_clipping(read)

			read = mask_low_quality_bases(read, seqQ_threshold)

			read = trim_3prime_N (read)

			# FUTURE: add option for low quality 3'trimming (for now there is only trimming of 3'Ns)

			subfamily [subidx] = read

		family_splitted [idx] = subfamily


	return family_splitted


def make_consensus_read (list_of_reads, method):
	'''
	STEPS 4 (single-strand consensus read) AND 5 (double-strand consensus read) OF PIPELINE
	Generates a consensus read from the given reads
	Requires: 
	(i) a list of read objects (pysam.AlignedSegment objects)
	(ii) specification on whether the consensus will be single-strand (method = "single_strand") or double (method = "double_strand")
	NOTE it also uses a verbose variable defined globally

	Returns: a single read object (pysam.AlignedSegment object), representing the consensus read of the input reads

	Sub-steps:
	1) Reconstruct alignment of the reads' sequence and base qualities, based on their CIGAR string. Please see "small" function "reconstruct_alignment" for more detailed information.
	2) Generate the consensus sequence and qualities. Please see "small" function "call_consensus" for more detailed information.
	3) Adjust fields (e.g. sequence, base qualities, position, CIGAR string, mapping quality) of the consensus read to the desired format (same format as pysam.AlignedSegment has)
	4) Create a read object (pysam.AlignedSegment object)

	NOTE: some parts should be improved in the future, such as calculating the consensus mapping quality (see below on comments of sub-step 3)
	The call_consensus() function might be also substituted for a more precise method when found.
	The "small functions" to assign ID, flag or tags to the consensus read might be also modified and customized in any desired way
	'''
	
	######   1. RECONSTRUCT ALIGNMENT   ######
	
	pos_reads, cigar_reads, seq_reads, qual_reads, mapQ_reads = extract_info_from_reads (list_of_reads)

	if verbose:
		print("\t\t * Reconstruct alignment in progress")

	aligned_seq, aligned_qual, minimum_pos = reconstruct_alignment (pos_reads, cigar_reads, seq_reads, qual_reads)

	
	#####    2. GENERATE CONSENSUS    ######

	if verbose:
		print("\t\t * Call consensus in progress") 

	aligned_cons_seq, aligned_cons_qual = call_consensus (aligned_seq, aligned_qual,
					      base_quality_shift,
					      max_base_quality,
					      error_rate_post_labeling,
					      error_rate_pre_labeling,
					      seqQ_threshold,
					      deletion_score,
					      no_insertion_score)
	

	#####    3. ADJUST CONSENSUS READ FIELDS to the desired format   ######

	if verbose:
		print("\t\t * Adjust consensus fields in progress") 

	cons_seq, cons_qual, cons_cigar, cons_pos = adjust_consensus_fields(aligned_cons_seq, aligned_cons_qual, minimum_pos)

	# Consensus mapping quality
	# it is calculated by average of all mapping qualities
	# therefore it should be improved
	# and this function might be joined to the previous function (adjust_consensus_fields), 
	# passing "mapQ_reads" as a new argument and returning a cons_mapQ value
	cons_mapQ = calculate_consensus_mapping_quality (mapQ_reads)

	# Consensus ID/query name
	# arbitrary reference read is given
	cons_id = get_consensus_id (list_of_reads[0], method)

	# Consensus flag
	# arbitrary reference read is given
	cons_flag = get_consensus_flag (list_of_reads[0], method)

	# Consensus chr/rname
	# arbitraty reference read is given
	# there has been a previous checkpoint where all reads from a family are checked to have the same rname
	cons_rname = list_of_reads[0].reference_id 

	# Tags - CHANGE - call function as optional - given argument. if false, just rx, mi - only if it takes a long time
	#cons_tags = (("MI", list_of_reads[0].get_tag("MI")), ("RX", list_of_reads[0].get_tag("RX")))
	cons_tags = add_tags(method, list_of_reads, aligned_seq, aligned_cons_seq, mapQ_reads)
	

	#####     4. CREATE READ OBJECT     #####
	
	cons_read = pysam.AlignedSegment()
	cons_read.query_name = cons_id
	cons_read.flag = cons_flag
	cons_read.reference_id = cons_rname
	cons_read.reference_start = cons_pos
	cons_read.mapping_quality = cons_mapQ
	cons_read.cigartuples = cons_cigar

	# mate info is not added here but in step function "fix_paired_end_fields()"

	cons_read.query_sequence = cons_seq
	cons_read.query_qualities = cons_qual
	cons_read.set_tags (cons_tags)
	
	return cons_read



def fix_paired_end_fields (paired_end_1, paired_end_2):
	'''
	Adds the mate/paired-end fields to the reads: RNEXT, PNEXT, TLEN

	Requires: two pysam.AlignedSegment objects, representing the two paired-end reads

	Returns: a list containing the same two pysam.AlignedSegment objects with the mate fieds added

	Assumes that the first read object corresponds to the left-most read (and vice-versa)
	'''

	rnames = [paired_end_1.reference_id,    paired_end_2.reference_id]
	pos    = [paired_end_1.reference_start, paired_end_2.reference_start]
	read_length =  paired_end_2.query_alignment_length

	# RNEXT (reference name - e.g. chr number)
	paired_end_1.next_reference_id = rnames[1]
	paired_end_2.next_reference_id = rnames[0]

	# PNEXT (starting position)
	paired_end_1.next_reference_start = pos[1]
	paired_end_2.next_reference_start = pos[0]

	# TLEN (template length)
	tlen = pos[1] + read_length - pos[0]
	paired_end_1.template_length = tlen
	paired_end_2.template_length = -tlen 


	return ([paired_end_1, paired_end_2])




#################################################           MAIN          ########################################################

def main():

	#########    PARSE ARGUMENTS     ########

	args = parse_args (sys.argv[1:])

	global input_filename
	input_filename = args.input_file

	global output_filename
	output_filename = args.output_file

	global verbose
	verbose = args.verbose

	global mapQ_threshold
	mapQ_threshold = args.min_map_quality
	
	global seqQ_threshold
	seqQ_threshold = args.min_base_quality

	global min_reads_for_consensus
	min_reads_for_consensus = args.min_reads

	global max_reads_for_consensus
	max_reads_for_consensus = args.max_reads
	
	global max_base_quality
	max_base_quality = args.max_base_quality

	global base_quality_shift
	base_quality_shift = args.base_quality_shift

	global error_rate_post_labeling
	error_rate_post_labeling = args.error_rate_post_labeling

	global error_rate_pre_labeling
	error_rate_pre_labeling = args.error_rate_pre_labeling

	global deletion_score
	deletion_score = args.deletion_score

	global no_insertion_score
	no_insertion_score = args.no_insertion_score


        ###########    OPEN FILES    #############

	###  Input file
	try:
		inbam = pysam.AlignmentFile(input_filename, 'rb')
		if verbose:
			print(input_filename, "has been read.")
	except:
		print ("ERROR: input file not found. \n Please specify file name of a valid .bam file. Include the format in the file name.")
		sys.exit(1)

	
	###  Output files
	
	if output_filename == None: 
		consensus_filename = "%s_cons.bam" % input_filename[:-4]	
	elif output_filename.endswith(".bam"):
		consensus_filename = output_filename
	else:
		print ("ERROR: output file is not specified in the right format. \n Please specify the file name of a valid .bam file. Include the format in the file name.")
		sys.exit(1)
 
	consensusbam = pysam.AlignmentFile(consensus_filename, 'wb', template = inbam)


	excluded_filename = "%s_filteredreads.bam" % consensus_filename[:-4]
	excludedbam = pysam.AlignmentFile(excluded_filename, 'wb', template = inbam)

	
	unprocessed_filename = "%s_filteredfamilies.bam" % consensus_filename[:-4]
	unprocessedbam = pysam.AlignmentFile(unprocessed_filename, 'wb', template = inbam)



        ############     DECLARE SOME VARIABLES     ###########

	passed_reads, excluded_reads, processed_families, excluded_families = 0, 0, 0, 0

	family, family_code, next_family = None, None, None

	global split_dict
	split_dict = {0: 'A1', 1: "B2", 2: "B1", 3: "A2"}



	############     READ BAM FILE     ###########

	for read in inbam:

                ###########      1. FILTER READS     ##########

		if pass_filters(read) == False:

			excluded_reads += 1
			excludedbam.write(read)

			continue 

		passed_reads += 1


                ###########      2. SELECT FAMILY    ##########

		family, family_code, next_family = add_read_to_family (read, family, family_code)

		if next_family == None:    # current family is not completely loaded in memory
					   # keep reading bam file 
			continue


		############      3. PREPROCESS      ############
				
		family_preprocessed = preprocess_family (family, family_code)

		if family_preprocessed == None: 	# not enough reads to form consensus

			excluded_families += 1

			for read_ in family:
				unprocessedbam.write(read_)

			family, family_code, next_family = load_next_family(next_family)

			continue
	
		processed_families += 1


		############    4. MAKE SINGLE-STRAND CONSENSUS SEQUENCE    ############

		single_consensus = [None, None, None, None]
		
		for idx, subfamily in enumerate(family_preprocessed):

			if verbose:
				print("Single-strand consensus for subfamily", family_code, split_dict[idx], "in progress")

			single_consensus[idx] = make_consensus_read (subfamily, method = "single_strand")
			


		############    5. MAKE DOUBLE CONSENSUS SEQUENCE    ############

		paired_end_1 = [single_consensus[0], single_consensus[1]] # consensus subfamilies A1, B2
		paired_end_2 = [single_consensus[2], single_consensus[3]] # consensus subfamilies B1, A2

		if verbose:
			print("Double-strand consensus for family", family_code, "in progress") 

		double_consensus = [make_consensus_read (paired_end_1, method = "double_strand"),
				    make_consensus_read (paired_end_2, method = "double_strand")]



		############    6. FIX PAIRED-END INFORMATION    ############ 
		
		fixed_consensus = fix_paired_end_fields (double_consensus[0], double_consensus[1])


		#########      WRITE IN BAM    ##########

		for consensus_read in fixed_consensus:
			consensusbam.write(consensus_read)

		if verbose:
			print("Consensus reads for family", family_code, "have been sucessfully written \n")

		####################################################################################


		#########      Update family    ##########

		family, family_code, next_family = load_next_family(next_family)



	#########      Repeat steps 3-6 for the last family    ########## 

	family_preprocessed = preprocess_family (family, family_code)
	if family_preprocessed == None:
		excluded_families += 1
		for read_ in family:
			unprocessedbam.write(read_)
	else:
		processed_families += 1
		single_consensus = [None, None, None, None]
		for idx, subfamily in enumerate(family_preprocessed):
			if verbose:
				print("Single-strand consensus for subfamily", family_code, split_dict[idx], "in progress") 
			single_consensus[idx] = make_consensus_read (subfamily, method = "single_strand")
		paired_end_1 = [single_consensus[0], single_consensus[1]]
		paired_end_2 = [single_consensus[2], single_consensus[3]]
		if verbose:
			print("Double-strand consensus for family", family_code, "in progress")
		double_consensus = [make_consensus_read (paired_end_1, method = "double_strand"), make_consensus_read (paired_end_2, method = "double_strand")]
		fixed_consensus = fix_paired_end_fields (double_consensus[0], double_consensus[1])
		for consensus_read in fixed_consensus:
			consensusbam.write(consensus_read)
		if verbose:
			print("Consensus reads for family", family_code, "have been sucessfully writen \n")

	if verbose:
		print ("\n Input file has been completely read \n")

	#####################################################################



       #########     CLOSE FILES    #########

	print("\n A total of %d reads (%.2f %%) passed the initial quality filters." % (passed_reads, passed_reads/(passed_reads+excluded_reads) * 100))
	print("\n A total of %d reads (%.2f %%) were filtered out." % (excluded_reads, excluded_reads/(passed_reads+excluded_reads) * 100))
	print("\n A total of %d families (%.2f %%) were successfully processed to generate a consensus read." % (processed_families, processed_families/(processed_families+excluded_families) * 100))
	print("\n A total of %d families (%.2f %%) were filtered out due to not enough reads to generate a consensus read." % (excluded_families, excluded_families/(processed_families+excluded_families) * 100))

	excludedbam   .close()
	unprocessedbam.close()
	consensusbam  .close()
	inbam         .close()




if __name__ == "__main__":
	main()
