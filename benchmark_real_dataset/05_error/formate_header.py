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
            m_size = re.search(r"(?P<size>size=\S+);", str(line))
            if m_size is not None:
                new_header.append(str(m_size.group('size')))
            m_count = re.search(r"(?P<count>count=\S+);", str(line))
            if m_count is not None:
                new_header.append(str(m_count.group('count')))
            m_sample = re.search(r"(?P<sample>merged_sample=\{.+\});",str(line))
            if m_sample is not None:
                new_header.append(str(m_sample.group('sample')).replace(" ", ""))
            print(">" + ";".join(new_header))
        else:
            print(line)
