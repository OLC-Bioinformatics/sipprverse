def test():
    mypileup = "/media/nas/akoziol/Pipeline_development/GeneSipperV2/baitTest/sequences/2015-SEQ-0795/Vibrio_MLST/recA/mpileup.out"
    with open(mypileup) as pileup:
        for line in pileup:
            if 'recA_97' in line:
                print line