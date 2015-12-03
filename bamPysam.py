
def test():
    import pysam
    import os
    # samfile = pysam.AlignmentFile("/media/nas/akoziol/Pipeline_development/GeneSipperV2/baitTest/sequences/2015-SEQ-0795/Vibrio_MLST/recA/2015-SEQ-0795_recA_bowtie2_sorted.bam", "rb" )
    samfile = pysam.AlignmentFile("/media/nas/akoziol/Pipeline_development/GeneSipperV2/baitTest/sequences/2015-SEQ-0795/Vibrio_MLST/recA/2015-SEQ-0795_recA_SMALT_sorted.bam", "rb" )

    # os.system('ls -larth /media/nas/akoziol/Pipeline_development/GeneSipperV2/baitTest/sequences/2015-SEQ-0795/Vibrio_MLST/recA/2015-SEQ-0795_recA_bowtie2_sorted.bam')

    # samfile = pysam.AlignmentFile("/media/nas/akoziol/Pipeline_development/GeneSipperV2/baitTest/sequences/2015-SEQ-0795/Vibrio_MLST/recA/2015-SEQ-0795_recA_bowtie2.sam", "r" )
    count = 0
    # for read in samfile.fetch('recA_97'):
    #     print read
        # count += 1
    # print count
    # for reference in samfile.references:
        # print reference
        # if reference == 'recA_97':
    for pileupcolumn in samfile.pileup('recA_97'):
        print pileupcolumn.pos, pileupcolumn.n
        for pileupread in pileupcolumn.pileups:
            print pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]
                # print ("\ncoverage at base %s = %s" %
                #     (pileupcolumn.pos, pileupcolumn.n))
        # for pileupread in pileupcolumn.pileups:
        #     if not pileupread.is_del and not pileupread.is_refskip:  # query position is None if is_del or is_refskip is set.
        #         print ('\tbase in read %s = %s' %
        #             (pileupread.alignment.query_name,
        #                  pileupread.alignment.query_sequence[pileupread.query_position]))

    # samfile.close()