from collections import defaultdict, Counter

def consensus_caller(input_reads, cutoff, tag, length_check):
    consensus_seq = ''
    if length_check is True:
        for read in input_reads[1:]:
            if len(read) != len(input_reads[0]):
                raise Exception((f"Read lengths for tag {tag} used for "
                                 f"calculating the SSCS are not uniform!!!"
                                 ))
    consensus_counter = Counter(input_reads).most_common(2)
    num_reads = len(input_reads)
    if len(consensus_counter) == 1:
        return consensus_counter[0][0]
    else:
        # In the order of T, C, G, A, N, Total
        nuc_key_dict = {0: 'T', 1: 'C', 2: 'G', 3: 'A', 4: 'N'}
        for i in range(len(input_reads[0])):  
            nuc_identity_list = [0, 0, 0, 0, 0, 0]
            # Count the types of nucleotides at a position in a read.
            # i is the nucleotide index within a read in groupedReadsList
            for j in range(num_reads):  
            # Do this for every read that comprises a tag family.
            # j is the read index within groupedReadsList
                try:
                    if input_reads[j][i] == 'T':
                        nuc_identity_list[0] += 1
                    elif input_reads[j][i] == 'C':
                        nuc_identity_list[1] += 1
                    elif input_reads[j][i] == 'G':
                        nuc_identity_list[2] += 1
                    elif input_reads[j][i] == 'A':
                        nuc_identity_list[3] += 1
                    elif input_reads[j][i] == 'N':
                        nuc_identity_list[4] += 1
                    else:
                        nuc_identity_list[4] += 1
                    nuc_identity_list[5] += 1
                except Exception:
                    break
            try:
                for j in [0, 1, 2, 3, 4]:
                    if (float(nuc_identity_list[j])
                            /float(nuc_identity_list[5])
                            ) >= cutoff:
                        consensus_seq += nuc_key_dict[j]
                        break
                    elif j == 4:
                        consensus_seq += 'N'
            except Exception:
                consensus_seq += 'N'
        return consensus_seq
