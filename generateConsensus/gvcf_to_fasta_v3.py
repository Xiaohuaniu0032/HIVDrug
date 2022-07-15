#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (C) 2021 Thermo Fisher Scientific Inc. All Rights Reserved
__version__ = '0.1.7'

import pysam
from optparse import OptionParser
import sys
import os
import json

# pysam compatibility
pysam_FastaFile = pysam.FastaFile if hasattr(pysam, 'FastaFile') else pysam.Fastafile
pysam_FastaFile_get_reference_length = lambda f_fasta, contig_name : (f_fasta.get_reference_length(contig_name) if hasattr(f_fasta, 'get_reference_length') else f_fasta.getReferenceLength(contig_name))  

def get_overlap(interval_1, interval_2):
    '''
    left-close, right-open intervals
    '''
    assert((len(interval_1) == len(interval_2) == 2) and (interval_1[0] <= interval_1[1]) and (interval_2[0] <= interval_2[1]))
    overlap_left = max(interval_1[0], interval_2[0])
    overlap_right = min(interval_1[1], interval_2[1])
    return (overlap_left, overlap_right) if (overlap_left <= overlap_right) else (None, None)

def is_overlap(interval_1, interval_2):
    '''
    left-close, right-open intervals
    '''
    my_overlap = get_overlap(interval_1, interval_2)
    if None in my_overlap:
        return False
        
    # Handle the case of my_overlap[0] == my_overlap[1]
    if my_overlap[0] == my_overlap[1]:
        # Case 1: [a, b) and [b, c) are not overlapping, where a<b<c
        # Case 2: [a, b) and [c, c) are overlapping and the overlap is [c, c), where a <= c <= b
        return set(my_overlap) in (set(interval_1), set(interval_2))

    return True

def is_too_far(allele_dict_1, allele_dict_2, max_pos_diff=256):
    if allele_dict_1['CHROM'] != allele_dict_2['CHROM']:
        return True
    interval_1 = (allele_dict_1['POS'], allele_dict_1['POS'] + len(allele_dict_1['REF']))
    interval_2 = (allele_dict_2['POS'], allele_dict_2['POS'] + len(allele_dict_2['REF']))
    if is_overlap(interval_1, interval_2):
        return False
    return max(interval_1[0], interval_2[0]) - min(interval_1[1], interval_2[1]) > max_pos_diff
    
    
def is_identical_variants(allele_dict_1, allele_dict_1_alt_idx, allele_dict_2, allele_dict_2_alt_idx, f_fasta, max_pos_diff=256):
    """
    Are the two allele identical (maybe with different alignment)?
    where allele_dict follows the VCF format (note that POS is 1-based), e.g.
    allele_dict = {'CHROM': 'chr1', 'POS': 65536, 'REF': 'TA', 'ALT': ['A'] 
    """

    # Trivial case 1: different chromosomes    
    if allele_dict_1['CHROM'] != allele_dict_2['CHROM']:
        return False
    
    # Flanking interval
    if allele_dict_1['POS'] < allele_dict_2['POS']:
        interval_left = (allele_dict_1['POS'], allele_dict_1['POS'] + len(allele_dict_1['REF']))
        interval_right = (allele_dict_2['POS'], allele_dict_2['POS'] + len(allele_dict_2['REF']))
        allele_dict_left = allele_dict_1
        allele_dict_right = allele_dict_2
    else:
        interval_left = (allele_dict_2['POS'], allele_dict_2['POS'] + len(allele_dict_2['REF']))
        interval_right = (allele_dict_1['POS'], allele_dict_1['POS'] + len(allele_dict_1['REF']))
        allele_dict_left = allele_dict_2
        allele_dict_right = allele_dict_1 
    
    flanking_start = interval_left[0]
    flanking_end = max(interval_left[1], interval_right[1])
    
    # Claim not identical if the two alleles span a too big interval
    if flanking_end - flanking_start > max_pos_diff:
        return False
    
    len_diff_1 = len(allele_dict_1['REF']) - len(allele_dict_1['ALT'][allele_dict_1_alt_idx]) 
    len_diff_2 = len(allele_dict_2['REF']) - len(allele_dict_2['ALT'][allele_dict_2_alt_idx]) 
    
    # Trivial case 2a: allele with different length change can't be identical    
    if len_diff_1 != len_diff_2:
        return False
    # Trivial case 2a: SNP with the same length
    elif len_diff_1 == 0 and len(allele_dict_1['REF']) == len(allele_dict_2['REF']) == 1:
        # The alleles are the same iff positions are the same.
        return allele_dict_1['POS'] == allele_dict_2['POS']

    # Trivial case 3: the same reference allele
    if allele_dict_1['POS'] == allele_dict_2['POS'] and allele_dict_1['REF'] == allele_dict_2['REF']:
        return allele_dict_1['ALT'][allele_dict_1_alt_idx] == allele_dict_2['ALT'][allele_dict_2_alt_idx]
    
    # Now get the flanking seq
    # No need to parse the fasta if the two alleles have overlapping.
    if is_overlap(interval_left, interval_right):
        num_suffix_from_right = flanking_end - interval_left[1]
        flanking_ref_seq = allele_dict_left['REF'] + allele_dict_right['REF'][-num_suffix_from_right:] if num_suffix_from_right else allele_dict_left['REF']
    else:
        # fasta fetch is 0-based
        flanking_ref_seq = f_fasta.fetch(allele_dict_1['CHROM'], flanking_start - 1, flanking_end - 1)        

    # flanking_alt_1
    num_left_paddings = allele_dict_1['POS'] - flanking_start
    num_right_paddings = flanking_end - (allele_dict_1['POS'] + len(allele_dict_1['REF']))
    if num_right_paddings:
        flanking_alt_1 = flanking_ref_seq[:num_left_paddings] + allele_dict_1['ALT'][allele_dict_1_alt_idx] + flanking_ref_seq[-num_right_paddings:]
    else:
        flanking_alt_1 = flanking_ref_seq[:num_left_paddings] + allele_dict_1['ALT'][allele_dict_1_alt_idx]
    # flanking_alt_2
    num_left_paddings = allele_dict_2['POS'] - flanking_start
    num_right_paddings = flanking_end - (allele_dict_2['POS'] + len(allele_dict_2['REF']))
    if num_right_paddings:
        flanking_alt_2 = flanking_ref_seq[:num_left_paddings] + allele_dict_2['ALT'][allele_dict_2_alt_idx] + flanking_ref_seq[-num_right_paddings:]
    else:
        flanking_alt_2 = flanking_ref_seq[:num_left_paddings] + allele_dict_2['ALT'][allele_dict_2_alt_idx]
    
    return flanking_alt_2 == flanking_alt_1

def parse_bed_record_basic(bed_line):
    strip_bed_line = bed_line.strip()
    split_bed_line = bed_line.split('\t')
    bed_dict = {
        'chrom': split_bed_line[0],
        'chromStart': int(split_bed_line[1]),
        'chromEnd': int(split_bed_line[2]),
        'name': split_bed_line[3],
        'RAWLINE': strip_bed_line
    }
    return bed_dict

def parse_vcf_record_heuristically(vcf_line):
    """
    Parse the VCF tags heursitcally. Didn't further split or convert INFO and FORMAT to the right type.
    """
    strip_vcf_line = vcf_line.strip()
    # Note that there is a bug that an SVB allele might contain empty REF or ALT. It could be wrong if I use  split_vcf_line = vcf_line.split().
    split_vcf_line = strip_vcf_line.split('\t')
    # INFO tags
    info_dict = dict([info_entry.split('=') if '=' in info_entry else (info_entry, None) for info_entry in split_vcf_line[7].split(';')])

    # QUAL
    try:
        qual = float(split_vcf_line[5])
    except ValueError:
        qual = split_vcf_line[5]

    # Create basic vcf dict
    vcf_dict = {
        'CHROM': split_vcf_line[0],
        'POS': int(split_vcf_line[1]),
        'ID': split_vcf_line[2],
        'REF': split_vcf_line[3],
        'ALT': split_vcf_line[4].split(','),
        'QUAL': qual,
        'FILTER': split_vcf_line[6],
        'INFO': info_dict,
        'FORMAT': [],
        'RAWLINE': strip_vcf_line,
    }

    # FORMAT tags
    if len(split_vcf_line) < 9:
        return vcf_dict
        
    format_dict_list = []
    format_key_list = split_vcf_line[8].split(':') if split_vcf_line[8] != '.' else []
    for sample_format_tag in split_vcf_line[9:]:
        split_sample_format_tag = sample_format_tag.split(':')
        if len(format_key_list) != len(split_sample_format_tag):
            format_dict_list.append({})
            continue
        
        format_dict = dict(zip(format_key_list, split_sample_format_tag))
        format_dict_list.append(format_dict)
        
    vcf_dict['FORMAT'] = format_dict_list
    
    return vcf_dict

def get_pos_end(vcf_dict):
    """
    return the 1-based, end position (included) of the vcf_dict
    """
    if is_coverage_gvcf_dict(vcf_dict):
        return int(vcf_dict['INFO']['END'])
    return vcf_dict['POS'] + len(vcf_dict['REF']) - 1

def ref_seq_from_coverage_gvcf_dict(gvcf_dict, min_dp, f_fasta):
    if not is_coverage_gvcf_dict(gvcf_dict):
        raise ValueError('The VCF line "%s" is not a coverage gVCF line.'%(gvcf_dict['RAWLINE']))
    my_dp = int(gvcf_dict['INFO']['DP'])
    # POS and END in gVCF is 1-based, left close, right close interval
    # pysam fasta fetch is 0-based, left close, right open interval
    chrom_start = int(gvcf_dict['POS']) - 1
    chrom_end = int(gvcf_dict['INFO']['END']) #- 1 + 1
    
    if my_dp < min_dp:
        return 'N'*(chrom_end - chrom_start)
    
    return f_fasta.fetch(gvcf_dict['CHROM'], chrom_start, chrom_end)
    
    
def calculate_vcf_gap_len(previous_vcf_dict, this_vcf_dict):
    if previous_vcf_dict['CHROM'] != this_vcf_dict['CHROM']:
        return None
    len_gap = this_vcf_dict['POS'] - get_pos_end(previous_vcf_dict) - 1
    return len_gap

def is_coverage_gvcf_dict(vcf_dict):
    return ''.join(vcf_dict['ALT']) == '.'


def do_variant_filtering_inplace(use_allele_idx_list, vcf_dict, min_non_hpindel_var_freq, min_hpindel_var_freq, min_dp):
    # Don't touch the special case of high QUAL LIA record
    if is_high_qual_lia_record(vcf_dict):
        return
    
    # Parse FDP
    fdp_str = vcf_dict['INFO'].get('FDP', vcf_dict['INFO'].get('DP', '.'))
    fdp = int(fdp_str) if fdp_str.isdigit() else 0
    # Filtering on coverage
    if fdp < min_dp:
        use_allele_idx_list[:] = []    
        vcf_dict['INFO']['FDP<MIN_DP'] = True
        return
    
    num_alt = len(vcf_dict['ALT'])
    fdvr_str = vcf_dict['INFO'].get('FDVR', ','.join('.'*num_alt))
    is_hp_indel = [fdvr == '0' for fdvr in fdvr_str.split(',')]
    af_str = vcf_dict['INFO'].get('AF', ','.join('0'*num_alt))
    af_list = [float(af) for af in af_str.split(',')]
    
    for allele_idx in tuple(use_allele_idx_list):
        if not allele_idx:
            continue
        alt_idx = allele_idx - 1
        
        my_min_var_freq_thresold = min_hpindel_var_freq if is_hp_indel[alt_idx] else min_non_hpindel_var_freq
        # Remove the alt_idx if its AF <  threshold
        if af_list[alt_idx] < my_min_var_freq_thresold:
            use_allele_idx_list.remove(allele_idx)

def is_high_qual_lia_record(vcf_dict):
    """
    Is the vcf_dict a high QUAL variant called by LIA?
    """
    return ('FDP' not in vcf_dict['INFO']) and vcf_dict['QUAL'] == 50.0 and (not is_coverage_gvcf_dict(vcf_dict))


def is_hp_diff(baseseq_0, baseseq_1):
    """
    Simple and effective (though not efficient) way to tell if the diff between two strings is just HP-INDEL
    """
    if len(baseseq_0) == len(baseseq_1):
        return False
    hp_compressed_list = [[], []]

    for seq_idx, baseseq in enumerate([baseseq_0, baseseq_1]):
        my_hp_compressed = hp_compressed_list[seq_idx]
        prev_nuc = None
        for nuc in baseseq:
            if nuc == prev_nuc:
                my_hp_compressed[-1][1] += 1
            else:
                my_hp_compressed.append([nuc, 1])
                prev_nuc = nuc
    if len(hp_compressed_list[0]) != len(hp_compressed_list[1]):
        return False
    
    num_nuc_diff = 0
    for nuc_order_idx in range(len(hp_compressed_list[0])):
        if hp_compressed_list[0][nuc_order_idx][0] != hp_compressed_list[1][nuc_order_idx][0]:
            return False
        if hp_compressed_list[0][nuc_order_idx][1] != hp_compressed_list[1][nuc_order_idx][1]:
            num_nuc_diff += 1
            if num_nuc_diff > 1:
                return False
    return True


def merge_tag_inplace(vcf_dict, info_key, source_idx, dest_idx, source_new_value=0, dtype=float):
    """
    merge the value of vcf_dict['INFO'][info_key][source_idx] to vcf_dict['INFO'][info_key][dest_idx]
    """
    my_tag = vcf_dict['INFO'].get(info_key, '.')
    if my_tag == '.':
        return
    
    split_tag = my_tag.split(',')
    if max(source_idx, dest_idx) >= len(split_tag) or min(source_idx, dest_idx) < 0:
        return
    
    my_tag_list = list(map(dtype, split_tag))
    my_tag_list[dest_idx] += my_tag_list[source_idx]
    my_tag_list[source_idx] = source_new_value
    new_tag = ','.join(map(str, my_tag_list))
    vcf_dict['INFO'][info_key] = new_tag
    
def heal_het_snp_complex_inplace(vcf_dict, f_fasta):
    """
    Special handling of het call consists of SNP + (SNP + small HP-error) with FDVR=10
    If detected, merge the AF of the error allele into SNP and modify GT
    """
    if is_coverage_gvcf_dict(vcf_dict):
        return
    
    gt_str = vcf_dict['FORMAT'][0].get('GT', None)
    if gt_str is None:
        return
    
    vcf_dict['FORMAT'][0].setdefault('ORIG_GT', gt_str)
    
    split_gt = gt_str.split('/')

    # Case of no heal: hom call, ref get called
    if len(set(split_gt)) < 2 or '0' in split_gt  or '.' in split_gt:      
        return
    
    # Get the allele idx in GT. Must be no '0' in GT
    alt_idx_pair = (int(split_gt[0]) - 1, int(split_gt[1]) - 1)
    gt_allele_0 = vcf_dict['ALT'][alt_idx_pair[0]]
    gt_allele_1 = vcf_dict['ALT'][alt_idx_pair[1]]

    # No heal if the len diff is too large
    if abs(len(gt_allele_0) - len(gt_allele_1)) not in (1, 2):      
        return
    
    # Check FDVR
    fdvr_list = vcf_dict['INFO'].get('FDVR', '.').split(',')
    if '.' in fdvr_list or len(fdvr_list) != len(vcf_dict['ALT']):
        return
    if fdvr_list[alt_idx_pair[0]] != '10' or fdvr_list[alt_idx_pair[1]] != '10':
        return
    
    # Check allele type :
    # The healing condition must be one snp and one other (i.e., complex) 
    gt_allele_type_pair = (determin_allele_type(vcf_dict['REF'], gt_allele_0), determin_allele_type(vcf_dict['REF'], gt_allele_1))
    if 'snp' not in gt_allele_type_pair:      
        return

    # We need one prefix/suffix padding to tell the hp diff
    extra_left_padding_pos0 = vcf_dict['POS'] - 2
    extra_right_padding_pos0 = vcf_dict['POS'] - 1 + len(vcf_dict['REF'])
    
    # Handle the pysam compatibility
    try:
        contig_len = pysam_FastaFile_get_reference_length(f_fasta, vcf_dict['CHROM'])
    except AttributeError:
        # In case I couldn't get the reference length.
        contig_len = extra_right_padding_pos0  
    
    extra_left_padding = f_fasta.fetch(vcf_dict['CHROM'], extra_left_padding_pos0, extra_left_padding_pos0 + 1) if extra_left_padding_pos0 >= 0 else ''
    extra_right_padding = f_fasta.fetch(vcf_dict['CHROM'], extra_right_padding_pos0, extra_right_padding_pos0 + 1) if extra_right_padding_pos0 < contig_len else ''
    padded_gt_allele_0 = extra_left_padding + gt_allele_0 + extra_right_padding
    padded_gt_allele_1 = extra_left_padding + gt_allele_1 + extra_right_padding

    # Check diff between the two alleles in GT    
    if not is_hp_diff(padded_gt_allele_0, padded_gt_allele_1):        
        return
        
    # All Pass: Now I start healing
    snp_idx = gt_allele_type_pair.index('snp')
    dest_allele_idx = alt_idx_pair[snp_idx]
    source_allele_idx = alt_idx_pair[not snp_idx]
    
    # Merge tags
    merge_tag_inplace(vcf_dict, 'AF', source_allele_idx, dest_allele_idx, 0.0, float)
    merge_tag_inplace(vcf_dict, 'AO', source_allele_idx, dest_allele_idx, 0, int)
    merge_tag_inplace(vcf_dict, 'FAO', source_allele_idx, dest_allele_idx, 0, int)

    # Update GT
    new_gt_str = '%s/%s' %(dest_allele_idx + 1, dest_allele_idx + 1)
    vcf_dict['FORMAT'][0]['GT'] = new_gt_str

    # Print the healing status in stdout
    healed_allele_str = vcf_dict_to_allele_str(vcf_dict, source_allele_idx  + 1)
    snp_allele_str = vcf_dict_to_allele_str(vcf_dict, dest_allele_idx + 1)
    print('The possibly erroneous allele %s is healed to become the SNP %s' %(healed_allele_str, snp_allele_str))

def determine_allele_idx_to_use(vcf_dict, major_allele_only):
    """
    Decide which allele to be written in the FASTA
    """
    gt_str = vcf_dict['FORMAT'][0]['GT']

    # REF allele
    if gt_str in ('0/0', './.', '0/.', './0'):
        return [0] 

    # HET case in Allele View VCF
    if gt_str in ('1/.', './1'): 
        return [1]

    gt_idx_list = sorted([int(my_gt_idx) for my_gt_idx in gt_str.split('/')])
    print("raw gt_idx_list is:",gt_idx_list)
    if len(gt_idx_list) != 2:
        raise ValueError('Unsupported GT: must be exactly two alleles.\n%s'% vcf_dict['RAWLINE'])
        
    # Simple case: HOM call
    if gt_idx_list[0] == gt_idx_list[1]:
        return [gt_idx_list[0]]         

    # Special logic: Always use the variant in LIA record with QUAL = 50
    if is_high_qual_lia_record(vcf_dict) and 0 in gt_idx_list:
        return [gt_idx_list[0]] if gt_idx_list[0] else [gt_idx_list[1]]

    # More complicated case: HET call
    allele_list = [vcf_dict['ALT'][gt_idx - 1] if gt_idx else vcf_dict['REF'] for gt_idx in gt_idx_list]
    print("allele_list is:",allele_list)
    is_same_length = len(allele_list[0]) == len(allele_list[1])
    print("is_same_length is:",is_same_length)
    # The Het alleles with the same lengths
    #if is_same_length and (not major_allele_only):
    #    return gt_idx_list
    
    if is_same_length and (not major_allele_only):
        # not major_allele_only
        num_alt = len(vcf_dict['ALT'])
        fdvr_str = vcf_dict['INFO'].get('FDVR', ','.join('.'*num_alt))
        is_hp_indel = [fdvr == '0' for fdvr in fdvr_str.split(',')]
        af_str = vcf_dict['INFO'].get('AF', ','.join('0'*num_alt))
        af_list = [float(af) for af in af_str.split(',')]

        
        fro = int(vcf_dict['INFO'].get('FRO', vcf_dict['INFO']['RO'])) # Flow Evaluator Reference allele observations
        fdp = int(vcf_dict['INFO'].get('FDP',vcf_dict['INFO']['DP']))  # Flow Evaluator read depth at the locus
        
        for allele_idx in tuple(gt_idx_list):
            if not allele_idx:
                # ref allele
                if fro == 0:
                    ref_freq = 0
                else:
                    ref_freq = float(fro)/fdp

                if ref_freq < min_non_hpindel_var_freq or ref_freq < min_hpindel_var_freq:
                    gt_idx_list.remove(allele_idx)

            else:
                # alt allele
                alt_idx = allele_idx - 1
                my_min_var_freq_thresold = min_hpindel_var_freq if is_hp_indel[alt_idx] else min_non_hpindel_var_freq
                # Remove the alt_idx if its AF <  threshold
                if af_list[alt_idx] < my_min_var_freq_thresold:
                    gt_idx_list.remove(allele_idx)

        return gt_idx_list
        

    # The Het alleles with the different lengths
    try:
        fro = int(vcf_dict['INFO'].get('FRO', vcf_dict['INFO']['RO']))
        fao_list = [int(fao) for fao in vcf_dict['INFO'].get('FAO', vcf_dict['INFO']['AO']).split(',')]
    except:
        print('WARNING: Bad FAO, AF, FRO, RO. The vcf line is skipped\n  %s'%vcf_dict['RAWLINE'])
        return [0]

    top_2_allele_cov = [fao_list[gt_idx - 1] if gt_idx else fro for gt_idx in gt_idx_list]
    best_allele_idx = gt_idx_list[0] if top_2_allele_cov[0] > top_2_allele_cov[1] else gt_idx_list[1]
    print("top_2_allele_cov is:",top_2_allele_cov)
    print("best_allele_idx is:",best_allele_idx)
    return [best_allele_idx]


def get_iupac_code(vcf_dict, use_allele_idx_list):
    if len(use_allele_idx_list) not in (0, 1, 2):
        raise ValueError('len(best_allele_idx_list) not in (0,1,2)\nbest_allele_idx_list=%s\n%s'% (use_allele_idx_list, vcf_dict['RAWLINE']))      
    
    # Put N for low coverage
    if vcf_dict['INFO'].get('FDP<MIN_DP', False):
        return 'N'*len(vcf_dict['REF'])
    
    # Empty use_allele_idx_list: possible all ALT got filtered out. Use REF instead.
    if not use_allele_idx_list:
        return vcf_dict['REF']
    
    # allele_idx = 0 for REF, i for vcf_dict['ALT'][i-1]
    allele_list = [vcf_dict['ALT'][allele_idx - 1] if allele_idx else vcf_dict['REF'] for allele_idx in use_allele_idx_list] 
    
    # Only one allele
    if len(use_allele_idx_list) == 1:
        return allele_list[0]
    
    
    if len(allele_list[0]) != len(allele_list[1]):
        raise ValueError('Two alleles in %s must has the same length' %allele_list)

    iupac_code_dict = {
        'AA': 'A',
        'TT': 'T',       
        'GG': 'G',
        'CC': 'C',
        'AT': 'W',
        'TA': 'W',
        'CG': 'S',
        'GC': 'S',
        'AC': 'M',
        'CA': 'M',
        'GT': 'K',
        'TG': 'K',
        'AG': 'R',
        'GA': 'R',
        'CT': 'Y',
        'TC': 'Y'
    }

    my_iupac_seq = ''.join([iupac_code_dict[''.join(zipped_bases)] for zipped_bases in zip(allele_list[0], allele_list[1])])
    return my_iupac_seq

def check_reference_genome(vcf_path, input_fasta_path, process_contig_name):
    fasta_basename_dict = {}
    reference_header_start = '##reference='
    contig_header_start = '##contig=<'
    all_contig_dict = {}
    use_fasta_path = None
    with open(vcf_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('#'):
                break
            if line.startswith('##reference='):
                # line should be, e.g., "##reference=/results/referenceLibrary/tmap-f3/hg19/hg19.fasta"
                ref_path = line[len(reference_header_start):]
                ref_basename = os.path.basename(ref_path)
                my_ref_path_list = fasta_basename_dict.setdefault(ref_basename, [])
                my_ref_path_list.append(ref_path)
            elif line.startswith(contig_header_start) and line.endswith('>'):
                # line should be "##contig=<ID=2019-nCoV,length=29903,assembly=Ion_AmpliSeq_SARS-CoV-2_reference>"
                one_contig_dict = dict([pair.split('=') for pair in line[len(contig_header_start):-1].split(',')])
                all_contig_dict[one_contig_dict['ID']] = one_contig_dict


    # Sanity Check: No reference fasta is found in the VCF header
    if not fasta_basename_dict:
        raise ValueError('No reference being found in the header of %s'%vcf_path)
    
    # Sanity Check: More than one reference is found in the VCF header: I don't know which one us actually used.
    if len(fasta_basename_dict.keys()) > 1:
        raise ValueError('Multiple reference genomes are found in the header of %s\n%s'%(vcf_path, ', '.join(fasta_basename_dict.keys())))
    

    # Sanity Check: process_contig_name must be in VCF header
    if process_contig_name.lower() != 'all' and (process_contig_name not in all_contig_dict):
        raise KeyError('process_contig_name=%s is not found in the contigs in the VCF header.'%process_contig_name)
    
    my_fasta_basename = fasta_basename_dict.keys()[0]
    
    # No input fasta is provided. Use the reference fasta in the header
    if input_fasta_path is None:
        ref_path_list = fasta_basename_dict[my_fasta_basename]
        for ref_path in ref_path_list:
            if os.path.exists(ref_path):
                use_fasta_path = ref_path
                break
        if use_fasta_path is None:
            raise IOError('The reference fasta %s in the VCF header does not exists.'%', '.join(ref_path_list))
    else:
        use_fasta_path = input_fasta_path
    
    # Check contig length
    f_fasta = pysam_FastaFile(use_fasta_path)
    for contig_name, contig_dict in all_contig_dict.items():
        if contig_name not in ('all', process_contig_name):
            continue
        fai_contig_len = pysam_FastaFile_get_reference_length(f_fasta, contig_name)
        vcf_contig_len = int(contig_dict['length'])
        if vcf_contig_len != fai_contig_len:
            raise ValueError('The contig length of %s does not match: %d in %s; %d in %s' %(contig_name, fai_contig_len, use_fasta_path, vcf_contig_len, vcf_path))
    f_fasta.close()
    
    # Sanity Check: The input fasta is different doesn't match the reference in the VCF header, which means that the VCF file is obtained by another reference gemome.        
    # I only display WARNING message instead of error out.
    if os.path.basename(use_fasta_path) not in fasta_basename_dict:
        print('WARNING: The input fasta %s is not found in the VCF header in %s, though the contig name and contig len matches.'%(use_fasta_path, vcf_path))

    return use_fasta_path

class FastaWriter:
    def __init__(self, fasta_w_path):
        self._f_fasta = open(fasta_w_path, 'w')
        self._fasta_batch = []
        self._fasta_batch_size = 0
        self._contig_name_list = []
        self._contig_stats_list = []        
        self.set_max_batch_size(65536)

    def __del__(self):
        if not self._f_fasta.closed:
            self.flush()            
            self._f_fasta.close()

    def __enter__(self):
        return self

    def __exit__(self, err_type, err_msg, err_traceback):
        if not self._f_fasta.closed:
            self.flush()            
            self._f_fasta.close()
        return False        

    def set_max_batch_size(self, max_batch_size):
        self._max_batch_size = max(0, int(max_batch_size))
        self.flush()
        return self._max_batch_size
    
    def flush(self):
        self._f_fasta.write(''.join(self._fasta_batch))
        self._fasta_batch = []
        self._fasta_batch_size = 0
        
    def write_contig(self, contig_name):
        if self._contig_name_list:
            if contig_name == self._contig_name_list[-1]:
                return
            elif contig_name in self._contig_name_list:
                raise ValueError('The contig_name "%s" has been written.' %(contig_name))
            header_start_str = '\n>'
        else:
            header_start_str = '>'
            
        self._contig_name_list.append(contig_name)
        self._contig_stats_list.append({'contig': contig_name, 'contig_len': 0, 'num_N': 0})            
        contig_header_line = header_start_str + contig_name + '\n'
        self.write_seq(contig_header_line, True)
        
    def write_seq(self, seq, is_header_line=False):
        if not self._contig_name_list:
            raise ValueError('At least one contig name must be written before writting the sequence.')
        self._fasta_batch.append(seq)
        self._fasta_batch_size += len(seq)
        if not is_header_line:
            self._contig_stats_list[-1]['contig_len'] += len(seq)
            self._contig_stats_list[-1]['num_N'] += seq.count('N')
        if self._fasta_batch_size > self._max_batch_size:
            self.flush()
        
    def get_contig_stats_tuple(self):
        return tuple(self._contig_stats_list)
    
def remove_duplicated_alleles_inplace(vcf_dict, prev_variants_list, f_fasta_in):
    """
    In case the same allele got reported twice by TVC and LIA, modify the "GT" in place if I found any duplicated variant being reported.
    """
    gt_str = vcf_dict['FORMAT'][0]['GT']
    vcf_dict['FORMAT'][0]['ORIG_GT'] = gt_str

    # REF allele
    if gt_str == '0/0' or '.' in gt_str:
        return
    
    gt_allele_idx_list = [int(allele_idx_str) if allele_idx_str.isdigit() else 0 for allele_idx_str in gt_str.split('/')]
    is_changed = False
    
    for prev_var_dict in prev_variants_list[::-1]:
       
        # break if too far (note that the for loop must iterate prev_variants_list from the back)
        if is_too_far(vcf_dict, prev_var_dict):
            break
        
        for gt_allele_idx_idx in range(len(gt_allele_idx_list)):
            gt_allele_idx = gt_allele_idx_list[gt_allele_idx_idx]

            # Skip ref allele
            if not gt_allele_idx:
                continue
            
            # check duplicated alleles
            if is_identical_variants(vcf_dict, gt_allele_idx - 1, prev_var_dict, 0, f_fasta_in):
                # Use ref if it is a dup allele
                gt_allele_idx_list[gt_allele_idx_idx] = 0
                is_changed = True
    
    if is_changed:
        vcf_dict['FORMAT'][0]['GT'] = '/'.join([str(my_idx) for my_idx in gt_allele_idx_list])


def dummy_vcf_dict():
    dummy_vcf_line = '\t'.join(['.', '0', '.', 'N', 'N', '.', '.', '.', '.', '.'])
    return parse_vcf_record_heuristically(dummy_vcf_line)

def extend_first_last_pseudo_gvcf_lines_in_place(sorted_gvcf_lines, sorted_region_bed_lines):
    """
    Add psudo-gVCF lines for writting "N" if gVCF does not full cover the designed region
    
    """
    # There is nothing I can help if no region BED file provided.
    if not sorted_region_bed_lines:
        return
    
    my_chrom = sorted_region_bed_lines[0]['chrom']
    if my_chrom != sorted_region_bed_lines[-1]['chrom']:
        raise ValueError('sorted_region_bed_lines must contain regions in one chromosome only.')

    min_chrom_start = sorted_region_bed_lines[0]['chromStart']
    max_chrom_end = max([region_dict['chromEnd'] for region_dict in sorted_region_bed_lines])
    
    # The case of no valid gVCF line: I will add one gVCF coverage track
    if not sorted_gvcf_lines:
        pseudo_gvcf_dict = dummy_vcf_dict()
        pseudo_gvcf_dict['CHROM'] = my_chrom
        pseudo_gvcf_dict['POS'] = min_chrom_start + 1
        pseudo_gvcf_dict['ALT'] = ['.']
        pseudo_gvcf_dict['INFO'] = {
            'DP': 0,
            'MAX_DP': 0,
            'MIN_DP': 0,
            'END': max_chrom_end #+ 1 - 1
        }
        sorted_gvcf_lines.append(pseudo_gvcf_dict)
        return

    # Check there is only one CHROM in sorted_gvcf_lines
    if sorted_gvcf_lines[0]['CHROM'] != sorted_gvcf_lines[-1]['CHROM']:
        raise ValueError('sorted_gvcf_lines must contain VCF lines in one chromosome only.')

    # Check consitency of chromosome
    first_gvcf_dict = sorted_gvcf_lines[0]    
    if first_gvcf_dict['CHROM'] != my_chrom:
        raise ValueError('Mismatched chromosomes in sorted_gvcf_lines and sorted_region_bed_lines')    
    
    # The case of the first vcf line comes after first region bed line (i.e., missing reads in the begining of the contig)
    if min_chrom_start + 1 < first_gvcf_dict['POS']:
        pseudo_gvcf_dict = dummy_vcf_dict()
        pseudo_gvcf_dict['CHROM'] = my_chrom
        pseudo_gvcf_dict['POS'] = min_chrom_start + 1
        pseudo_gvcf_dict['ALT'] = ['.']
        pseudo_gvcf_dict['INFO'] = {
            'DP': 0,
            'MAX_DP': 0,
            'MIN_DP': 0,
            'END': first_gvcf_dict['POS'] - 1
        }
        sorted_gvcf_lines.insert(0, pseudo_gvcf_dict)
    
    # The case of last vcf line comes after first region bed line
    last_gvcf_dict = sorted_gvcf_lines[-1]
    last_gvcf_pos_end = get_pos_end(last_gvcf_dict)
    if  last_gvcf_pos_end < max_chrom_end:
        pseudo_gvcf_dict = dummy_vcf_dict()
        pseudo_gvcf_dict['CHROM'] = my_chrom
        pseudo_gvcf_dict['POS'] = last_gvcf_pos_end + 1
        pseudo_gvcf_dict['ALT'] = ['.']
        pseudo_gvcf_dict['INFO'] = {
            'DP': 0,
            'MAX_DP': 0,
            'MIN_DP': 0,
            'END': max_chrom_end
        }
        sorted_gvcf_lines.append(pseudo_gvcf_dict)


def vcf_dict_to_allele_str(vcf_dict, allele_idx):
    """
    Generate a string to represent the ALT allele
    """
    alt_idx = allele_idx - 1
    if alt_idx < 0:
        return '%s:%s_%s' %(vcf_dict['CHROM'], vcf_dict['POS'], vcf_dict['REF'])

    return '%s:%s_%s_%s' %(vcf_dict['CHROM'], vcf_dict['POS'], vcf_dict['REF'], vcf_dict['ALT'][alt_idx])

def calculate_num_vcf_paddings(ref, alt):
    """
    Determine the prefix, suffix paddings for the ref/alt in VCF format (i.e., empty allele is not allowed)
    """
    num_right_padding = 0 
    num_left_padding = 0
    min_len = min(len(ref), len(alt))
    
    # Trivial case
    if min_len == 1:
        return num_left_padding, num_right_padding

    # Assume left alignment. So I remove suffix padding first.
    while num_right_padding + 1 < min_len:
        if ref[-(num_right_padding + 1)] == alt[-(num_right_padding + 1)]:
            num_right_padding += 1
        else:
            break

    # Remove suffix padding
    if num_right_padding > 0:
        min_len -= num_right_padding
    
    # Then I remove the prefix padding.
    while num_left_padding + 1 < min_len:
        if ref[num_left_padding] == alt[num_left_padding]:
            num_left_padding += 1
        else:
            break
    
    return num_left_padding, num_right_padding

def determin_allele_type(ref, alt):
    """
    Determine the allele type (similar to the logic of generate_variant_table.py)
    """
    # SNP or MNP
    if len(ref) == len(alt):
        num_diff = len([1 for idx in range(len(ref)) if ref[idx] != alt[idx]])
        return 'snp' if num_diff == 1 else 'other'

    # Calculate padding (in VCF format)
    num_left_padding, num_right_padding = calculate_num_vcf_paddings(ref, alt)

    # Strip paddings
    if num_right_padding:
        ref = ref[num_left_padding:-num_right_padding]
        alt = alt[num_left_padding:-num_right_padding]
    elif num_left_padding:
        ref = ref[num_left_padding:]
        alt = alt[num_left_padding:]

    if alt.startswith(ref):
        return 'ins'
    elif ref.startswith(alt):
        return 'del'
    return 'other'
    
    
def update_summary_dict_inplace(summary_dict, vcf_dict, use_allele_idx_list):
    """
    Update the content of summary_dict for the processing status of vcf_dict
    """
    orig_gt = vcf_dict['FORMAT'][0].get('ORIG_GT', vcf_dict['FORMAT'][0].get('GT', '.'))
    if orig_gt in ('.', '0/0', './.'):
        return
        
    split_orig_gt = orig_gt.split('/')
    zygosity_key = 'homo' if len(set(split_orig_gt)) == 1 else 'het'

    for a_idx in set(split_orig_gt):
        if a_idx in ('0', '.'):
            continue
        a_idx = int(a_idx)
        my_allele_str = vcf_dict_to_allele_str(vcf_dict, a_idx)
        if a_idx in use_allele_idx_list:
            my_writte_key = 'variants_written'
            allele_type = determin_allele_type(vcf_dict['REF'], vcf_dict['ALT'][a_idx-1])
            if allele_type in ('ins', 'del'):
                allele_type = 'indel'
            my_variant_key = zygosity_key + '_' + allele_type + 's' if allele_type != 'other' else allele_type 
            summary_dict['variants_by_chromosome'][0][my_variant_key] += 1
            summary_dict['variants_by_chromosome'][0]['variants'] += 1
        else:
            my_writte_key = 'variants_not_written'
        summary_dict[my_writte_key].append(my_allele_str)

def gvcf_to_fasta_main(gvcf_path, input_fasta_path, output_fasta_path, process_contig_name, **kwargs):
    """
    Main function to convert gVCF to FASTA
    """
    # Parse options
    input_fasta_file = kwargs.get('input_fasta_file', None)
    region_bed_file = kwargs.get('region_bed_file', None)
    alias_contig = kwargs.get('alias_contig', process_contig_name)
    min_dp = kwargs.get('min_dp', 20)
    major_allele_only = kwargs.get('major_allele_only', 1)
    min_non_hpindel_var_freq = kwargs.get('min_non_hpindel_var_freq', 0.6)
    min_hpindel_var_freq = kwargs.get('min_hpindel_var_freq', 0.5)
    summary_json_path = kwargs.get('summary_json_path', os.path.join(os.path.dirname(output_fasta_path), 'summary.json') )
    
    
    use_fasta_path = check_reference_genome(gvcf_path, input_fasta_file, process_contig_name)

    if input_fasta_path is None:
        print('No reference FASTA provided. Use the reference FASTA specified in the VCF header: %s'%use_fasta_path)

    # alias_contig_dict is the map from original contig name to the contig name in the output FASTA
    # Although I support single contig only, I reserve the multi-contig for future expansion.
    alias_contig_dict = {process_contig_name: alias_contig}
    
    # Initialize summary_dict
    variants_by_chromosome_dict = {
        'chromosome'    : process_contig_name,
        'variants'      : 0,
        'het_snps'      : 0,
        'homo_snps'     : 0,
        'het_indels'    : 0,
        'homo_indels'   : 0,
        'other'         : 0        
    }
    meta_dict = {
        'version': __version__,
        'gvcf_path': os.path.abspath(gvcf_path),
        'process_contig': process_contig_name,
        'alias_contig': alias_contig,
        'region_bed_file': '' if region_bed_file is None else os.path.abspath(region_bed_file),
        'reference_fasta_path': os.path.abspath(use_fasta_path),
        'output_fasta_path': os.path.abspath(output_fasta_path)
    }
    param_dict = {
        'min_dp': min_dp,
        'major_allele_only': major_allele_only,
        'min_non_hpindel_var_freq': min_non_hpindel_var_freq,
        'min_hpindel_var_freq': min_hpindel_var_freq            
    }
    summary_dict = {
        'variants_written': [],
        'variants_not_written': [],
        'variants_by_chromosome': [variants_by_chromosome_dict],
        'meta': meta_dict,
        'parameters': param_dict
    }
    
    # Parse gVCF lines for the contig that I want to process
    with open(gvcf_path, 'r') as f_gvcf:
        gvcf_lines = [parse_vcf_record_heuristically(line) for line in f_gvcf if line.startswith(process_contig_name + '\t')]

    # Sort gvcf_lines by POS
    gvcf_lines.sort(key=lambda my_dict: my_dict['POS'])

    
    # Parse region BED file
    region_bed_lines = []
    if region_bed_file is not None:
        with open(region_bed_file, 'r') as f_bed:
            region_bed_lines = [parse_bed_record_basic(line) for line in f_bed if line.startswith(process_contig_name + '\t')]
        # Sort region_bed_lines by chromStart
        region_bed_lines.sort(key=lambda my_dict: my_dict['chromStart'])

    # Extend pseudo-gVCF lines in case no read covered in the first/last amplicon
    extend_first_last_pseudo_gvcf_lines_in_place(gvcf_lines, region_bed_lines)

    
    with FastaWriter(output_fasta_path) as f_fasta_out:
        # Old pysam doesn't support with statenent for pysam.Fastafile
        f_fasta_in = pysam_FastaFile(use_fasta_path)

        # prev_most_right_gvcf_dict is the gVCF line that covered the most right position previously
        prev_most_right_gvcf_dict = dummy_vcf_dict()
        # prev_variants_list is the list of dict that stored the variants being used for FASTA
        prev_variants_list = []
        
        for gvcf_dict in gvcf_lines:
            # Not the contig of interests, though shouldn't happen in this version of code
            if process_contig_name != gvcf_dict['CHROM']:
                prev_most_right_gvcf_dict = gvcf_dict
                continue

            # Write contig anyway (f_fasta_out knows how to deal with old or new contig)
            f_fasta_out.write_contig(alias_contig_dict[gvcf_dict['CHROM']])
            
            # Handle overlapping VCF record
            prev_most_right_end_pos = get_pos_end(prev_most_right_gvcf_dict)
            # Is the cuurent gVCF line covered (overlapped) with any gVCF line that has be processed previously?
            print(gvcf_dict)
            print("prev_most_right_gvcf_dict")
            print(prev_most_right_gvcf_dict)
            print("\n")
            is_covered_by_prev_right_most = gvcf_dict['CHROM'] == prev_most_right_gvcf_dict['CHROM'] and prev_most_right_end_pos >= gvcf_dict['POS']
            print(is_covered_by_prev_right_most)
            # If yes, then trim the position of the current gVCF line. Otherwise, the same position could be written multiple times to the FASTA.
            if is_covered_by_prev_right_most:
                # Trim the current gVCF coverage line if it is not fully covered by prev_most_right_end_pos
                if is_coverage_gvcf_dict(gvcf_dict) and prev_most_right_end_pos < get_pos_end(gvcf_dict):
                    trimmed_gvcf_dict = dummy_vcf_dict()
                    trimmed_gvcf_dict['CHROM'] = gvcf_dict['CHROM']
                    trimmed_gvcf_dict['POS'] = prev_most_right_end_pos + 1
                    trimmed_gvcf_dict['REF'] = 'N'
                    trimmed_gvcf_dict['ALT'] = ['.']
                    trimmed_gvcf_dict['INFO'] = {'DP': gvcf_dict['INFO']['DP'], 'END': gvcf_dict['INFO']['END']}
   
                    my_seq = ref_seq_from_coverage_gvcf_dict(trimmed_gvcf_dict, min_dp, f_fasta_in)
                    f_fasta_out.write_seq(my_seq)
                # Drop out the variant VCF line if overlap or it is fully covered by prev_most_right_end_pos
                else:
                    update_summary_dict_inplace(summary_dict, gvcf_dict, [])
                    print('WARNING: Overlapping VCF lines! The second line is discarded!\n  %s\n  %s' %(prev_most_right_gvcf_dict['RAWLINE'], gvcf_dict['RAWLINE']))

                # Do update here because I will continue soon
                if prev_most_right_end_pos < get_pos_end(gvcf_dict):
                    prev_most_right_end_pos = gvcf_dict
                continue
            
            # Calculate coverage gap
            cov_gap = calculate_vcf_gap_len(prev_most_right_gvcf_dict, gvcf_dict)
            if cov_gap > 0:
                f_fasta_out.write_seq('N'*cov_gap)
                
            # Get sequence from gVCF coverage line
            if is_coverage_gvcf_dict(gvcf_dict):
                # Put REF sequence or N's if no variant discovered
                my_seq = ref_seq_from_coverage_gvcf_dict(gvcf_dict, min_dp, f_fasta_in)

            # Get sequence from variant calling line
            else:
                # Revert a called ALT allele to REF if the ALT has been reported previousely.
                # Example: A right-aliged INS called by TVC but the left-aligned INS called by LIA
                remove_duplicated_alleles_inplace(gvcf_dict, prev_variants_list, f_fasta_in)                

                # Heal snp/complex
                heal_het_snp_complex_inplace(gvcf_dict, f_fasta_in)
                
                # Determine the allele idx to be used for the FASTA sequence
                use_allele_idx_list = determine_allele_idx_to_use(gvcf_dict, major_allele_only)
                print("use allele idx list is:",use_allele_idx_list)

                # Do allele filtering inplace
                do_variant_filtering_inplace(use_allele_idx_list, gvcf_dict, min_non_hpindel_var_freq, min_hpindel_var_freq, min_dp)

                # Get the varaint sequence (possibly with IUPAC code WSKMYR)
                my_seq = get_iupac_code(gvcf_dict, use_allele_idx_list)

                # Keep the variants that were output in FASTA. They will be used for remove_duplicated_alleles_inplace.
                for allele_idx in use_allele_idx_list:
                    # Skip REF allele
                    if not allele_idx:
                        continue

                    tmp_variant_dict = {
                        'CHROM': gvcf_dict['CHROM'],
                        'POS': gvcf_dict['POS'],
                        'REF': gvcf_dict['REF'],
                        'ALT': [gvcf_dict['ALT'][allele_idx - 1]],
                    }
                    prev_variants_list.append(tmp_variant_dict)

                # Update alt written status dict
                update_summary_dict_inplace(summary_dict, gvcf_dict, use_allele_idx_list)

            # Write the sequence to output FASTA
            f_fasta_out.write_seq(my_seq)
            
            # Do updates
            if not is_covered_by_prev_right_most:
                prev_most_right_gvcf_dict = gvcf_dict

        f_fasta_in.close()
        
        # Get contig_stats_tuple
        summary_dict['output_contig_stats'] = list(f_fasta_out.get_contig_stats_tuple())
    
    # Write summary.json    
    json.dump(summary_dict, open(summary_json_path, 'w'), indent=4, sort_keys=1)

    
if __name__ == '__main__':
    
    # Option Parser
    parser = OptionParser('Generate consensus FASTA from variant calling results.')
    parser.add_option('-v', '--gvcf-file',                help='(Required) Genome VCF (gVCF) file obtained from variant caller', dest='gvcf_file')
    parser.add_option('-f', '--input-fasta-file',         help='The reference genome FASTA file used for variant calling [default=the reference genome specified in the gVCF header]', dest='input_fasta_file')
    parser.add_option('-r', '--region-bed-file',          help='(Optional) Region BED file used to improve the reporting the "N" in the first and last regions.', dest='region_bed_file')
    parser.add_option('-o', '--output-fasta-file',        help='(Required) Output gemome FASTA file', dest='output_fasta_file')
    parser.add_option('-d', '--min-dp',                   help='Minimum DP required for trusting the original reference [default=20]', dest='min_dp')
    parser.add_option('-c', '--process-contig',           help='(Required) Process the contig in input_fasta_file only', dest='process_contig')
    parser.add_option('-a', '--alias-contig',             help='(Optional) The contig name in the output FASTA [default is the same as process_contig]', dest='alias_contig')
    parser.add_option('-m', '--major-allele-only',        help='(Optional) (0 or 1) Use the allele of major allele frequency only [default=1]', dest='major_allele_only')
    parser.add_option('-p', '--min-hpindel-var-freq',     help='(Optional) The minimum variant frequency of a HP-INDEL to put it in the FASTA. [default=0.6]', dest='min_hpindel_var_freq')
    parser.add_option('-n', '--min-non-hpindel-var-freq', help='(Optional) The minimum variant frequency of a variant that is not a HP-INDEL to put it in the FASTA. [default=0.5]', dest='min_non_hpindel_var_freq')


    if len(sys.argv) == 1 or '--help' in sys.argv:
        parser.print_help()
        sys.exit(1)
    (options, args) = parser.parse_args()
    
    
    # Parse gvcf_file gvcf
    if options.gvcf_file is not None:
        gvcf_file = options.gvcf_file
    else:
        raise IOError('Genome VCF file is not specified by "-v" or "--gvcf-file"')

    # Parse input_fasta_file
    input_fasta_file = options.input_fasta_file
    
    # Parse region_bed_file
    region_bed_file = options.region_bed_file
    
    # Parse output_fasta_file    
    if options.output_fasta_file is not None:
        output_fasta_file = options.output_fasta_file
    else:
        raise IOError('Output FASTA file is not specified by "-o" or "--output-fasta-file"')
    
    # Parse min_dp
    min_dp = int(options.min_dp) if options.min_dp is not None else 20
    
    # Parse process_contig
    if options.process_contig is not None:
        process_contig = options.process_contig
    else:
        raise IOError('The contig to be processed must be specified via "-c" or "--process-contig" option')

    # Parse alias_contig
    alias_contig = options.alias_contig if options.alias_contig is not None else options.process_contig
    
    # Parse major_allele_only
    major_allele_only = options.major_allele_only if options.major_allele_only is not None else '1'
    major_allele_only = int(major_allele_only)
    if major_allele_only not in (0, 1):
        raise ValueError('The value of "-m" or "--major-allele-only" must be either 0 or 1.')
    
    # Parse min_hpindel_var_freq
    min_hpindel_var_freq = options.min_hpindel_var_freq if options.min_hpindel_var_freq is not None else 0.6
    try:
        min_hpindel_var_freq = float(min_hpindel_var_freq)
    except ValueError:
        raise ValueError('The value of "-p" or "--min-hpindel-var-freq" must be a floating point number.')
    if not (0.0 <= min_hpindel_var_freq <= 1.0):
        raise ValueError('The value of "-p" or "--min-hpindel-var-freq" must be between 0 and 1.')        

    # Parse min_non_hpindel_var_freq
    min_non_hpindel_var_freq = options.min_non_hpindel_var_freq if options.min_non_hpindel_var_freq is not None else 0.5
    try:
        min_non_hpindel_var_freq = float(min_non_hpindel_var_freq)
    except ValueError:
        raise ValueError('The value of "-n" or "--min-non-hpindel-var-freq" must be a floating point number.')
    if not (0.0 <= min_non_hpindel_var_freq <= 1.0):
        raise ValueError('The value of "-n" or "--min-non-hpindel-var-freq" must be between 0 and 1.')        
    
    # Path to summary.json
    summary_json_path = os.path.join(os.path.dirname(output_fasta_file), 'summary.json')
    
    # Main function
    gvcf_to_fasta_main(gvcf_file, input_fasta_file, output_fasta_file, process_contig,
                       alias_contig=alias_contig,
                       min_dp=min_dp,
                       min_non_hpindel_var_freq=min_non_hpindel_var_freq,
                       min_hpindel_var_freq=min_hpindel_var_freq,
                       major_allele_only=major_allele_only,
                       input_fasta_file=input_fasta_file,
                       region_bed_file=region_bed_file,
                       summary_json_path=summary_json_path
    )

    print('The FASTA file is successfully written to "%s".' %output_fasta_file)
    print('Summary is written to "%s"'%summary_json_path)
