from collections import namedtuple
import traceback
from pyfasta import Fasta
from fasta import align, smatrix
import progressbar


library_path = 'lib.fasta'
query_path = 'Query.txt'
library_process_limit =  None #100

AlignmentResult = namedtuple('Result', ['title', 'score'])

'''
# Creates alignment outline
'''
def outline_alignment_for(aligned1, aligned2):
    alignment = ''

    for i in range(len(aligned1)):
        if aligned1[i] == aligned2[i]:
            alignment += '|'
        else:
            if aligned1[i] == '-' or aligned2[i] == '-':
                alignment += ' '
            else:
                alignment += '.'

    return alignment


'''
# Note that we should explicitly convert pyfasta record to string, as, otherwise,
# it will never be read into the memory, what will result in poor performance
# due to disk access.
'''
def process_query():
    print('Reading sequence library and query sequence')
    library = Fasta(library_path)
    queries = Fasta(query_path)
    query_sequence = str(queries["Rattus"])

    print('Processing')
    progress = progressbar.ProgressBar()
    results = []

    for record in progress(list(library.keys())[:library_process_limit]):
        try:
            library_sequence = str(library[record])
            aligned1, aligned2, score = align(library_sequence, query_sequence)
            results.append(AlignmentResult(title=record, score=score))
        except Exception as e:
            print(e.__class__, e)
            traceback.print_tb(e.__traceback__)
            print("Failed on "+ record)
            exit()

    etalone_score = sum([ smatrix[(x, x)] for x in query_sequence ])

    print("Done")
    print("Etalone score is %d" % etalone_score)
    print("Got %d results, here are top-30 among them:" % len(results))
    print("Score  | Match  | Record")

    for sequence in sorted(results, key=lambda x: x.score, reverse=True)[:30]:
        match = (sequence.score / etalone_score) * 100.0
        print("%6d | %5.3f%% | %s" % (sequence.score, match, sequence.title))

def test_simple():
    #str1 = "TACCGA"
    #str2 = "ACTGAC"

    str1 = "HEAGAWGHEE"
    str2 = "PAWHEAE"

    aligned1, aligned2, score = align(db_seq=str1, query_seq=str2)
    alignment = outline_alignment_for(aligned1, aligned2)
    print(aligned1)
    print(alignment)
    print(aligned2)
    print(score)

def test_rattus():
    library = Fasta(library_path)
    queries = Fasta(query_path)

    #query_sequence = str(queries["Rattus"])
    #query_sequence = str(library["sp|P06757|ADH1_RAT Alcohol dehydrogenase 1 OS=Rattus norvegicus GN=Adh1 PE=1 SV=3"])
    query_sequence = str(library["tr|G0Z9N0|G0Z9N0_MALDO Alcohol dehydrogenase (Fragment) OS=Malus domestica GN=Adh1-1 PE=4 SV=1"])
    library_sequence = str(library["tr|Q8K571|Q8K571_RAT ADH-like protein OS=Rattus norvegicus GN=Adh1 PE=2 SV=1"])

    aligned1, aligned2, score = align(library_sequence, query_sequence)
    outline = outline_alignment_for(aligned1, aligned2)

    with open('rattus_test.txt', 'a') as fl:
        print(query_sequence, file=fl)
        print(library_sequence, file=fl)
        print(score, file=fl)
        print(aligned1, file=fl)
        print(outline,  file=fl)
        print(aligned2, file=fl)


def main():
    #process_query()
    #test_simple()
    test_rattus()

if __name__ == "__main__":
    main()

