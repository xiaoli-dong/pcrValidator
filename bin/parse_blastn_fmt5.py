from Bio.Blast import NCBIXML
from Bio.Blast import Record
from Bio.Blast.Record import HSP

# https://uoftcoders.github.io/studyGroup/lessons/python/biopython/lesson/

result_handle = open("template.out.5.txt", "r")
blast_records = NCBIXML.parse(result_handle)
E_VALUE_THRESH = 10

seen = []

for blast_record in blast_records:
    # print("num_alignments:", blast_record.num_alignments)
    # query id
    # print(blast_record.query)
    # query lenght
    # print(blast_record.query_letters)
    qlen = blast_record.query_letters

    for alignment in blast_record.alignments:

        for hsp in alignment.hsps:
            align_len = hsp.align_length
            if hsp.expect < E_VALUE_THRESH and (align_len / qlen > 0.90):
                # print("****Alignment****")
                hit_def = alignment.hit_def

                target_id = hit_def.split(" ", 1)[0]
                id = f"{target_id}_{hsp.sbjct_start}_{hsp.sbjct_end}_{hsp.strand[1]} qlen={qlen} qstart={hsp.query_start} qend={hsp.query_end}"

                # print("sequence:", alignment.title)
                # print("id:", alignment.hit_id)
                if id not in seen:
                    seen.append(id)
                    print(f">{id}")
                    # print("length:", alignment.length)
                    # print("e value:", hsp.expect)
                    # print(hsp.query[0:75] + "...")
                    # print(hsp.match[0:75] + "...")
                    # print(hsp.sbjct[0:75] + "...")
                    # print(hsp.query + "...")
                    # print(hsp.match + "...")
                    print(hsp.sbjct)
                    # print("align length:", hsp.align_length)
                    # print("strand:", hsp.strand)
                    # print("evalue:", hsp.expect)
                    # print("identities:", hsp.identities)
