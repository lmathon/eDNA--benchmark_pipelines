#!/usr/bin/env python3

import sys, re

if __name__ == "__main__":
    for line in sys.stdin:
        line = line.strip()
        if line.startswith('>'):
            new_header = []
            m_seq = re.search(r"^>(?P<seq_id>\S+).+", str(line))
            if m_seq is not None:
                new_header.append(str(m_seq.group('seq_id')))
            m_family = re.search(r"(?P<family>family_name=\S+);", str(line))
            if m_family is not None:
                new_header.append(str(m_family.group('family')))
            m_genus = re.search(r"(?P<genus>genus_name=\S+);", str(line))
            if m_genus is not None:
                new_header.append(str(m_genus.group('genus')))
            m_species = re.search(r"(?P<species>species_name=\D+);",str(line))
            if m_species is not None:
                new_header.append(str(m_species.group('species')))
            #m_rank = re.search(r"(?P<rank>rank=\S+);",str(line))
            #if m_rank is not None:
            #    new_header.append(str(m_rank.group('rank')))
            print(">" + ";".join(new_header)+";")
        else:
            print(line)
