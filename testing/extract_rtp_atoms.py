#!/usr/bin/python3

rtp_file = "/home/joaov/python-mmpbsa/mmpbsa/testing/ffG54a7pHt.rtp"
lines = [] 
found_section=False
exclusion_list = ('[ bonds ]','[ angles ]','[ impropers ]','[ dihedrals ]','[ exclusions ]',
                ' [ bonds ]',' [ angles ]',' [ impropers ]',' [ dihedrals ]',' [ exclusions ]')
aa_list = []
with open(rtp_file) as rtp:
    next(rtp)
    next(rtp)
    next(rtp)
    for line in rtp:
        if line.startswith('[ ') and not line.startswith(exclusion_list):
            aa_list.append(line)
            found_section = True
            cline = str(line).rstrip('\n')
            lines.append(cline)
            continue
        if found_section:
            if line.startswith(exclusion_list):
                found_section = False
            elif '[ atoms ]' in line:
                continue
            else:
                s_line = str(line).rstrip('\n')
                lines.append(s_line)

log_rtp = "/home/joaov/python-mmpbsa/mmpbsa/testing/log_rtp.rtp"
with open(log_rtp, "w") as log:
    log.writelines("%s\n" % place for place in lines)