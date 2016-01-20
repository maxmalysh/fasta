from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
from pyfasta import Fasta
import traceback
import progressbar

from fasta import align, smatrix, get_performance_timer
from utils import outline_alignment_for, SynchronizedTimer

library_path = 'lib.fasta'
query_path = 'Query.txt'
library_process_limit = None

AlignmentResult = namedtuple('Result', ['title', 'score'])


class AlignmentTask:
    def __init__(self, record, future):
        self.record = record
        self.future = future
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
    progress = progressbar.ProgressBar(max_value=len(library.keys()))
    cpu_count = multiprocessing.cpu_count()
    executor = ThreadPoolExecutor(max_workers=cpu_count)

    tasks = []
    for record in list(library.keys())[:library_process_limit]:
        library_sequence = str(library[record])
        future = executor.submit(align, library_sequence, query_sequence)
        tasks.append(AlignmentTask(record, future))

    results = []
    for i in range(len(tasks)):
        _, _, score = tasks[i].future.result()
        results.append(AlignmentResult(title=tasks[i].record, score=score))
        progress.update(i)

    etalone_score = sum([ smatrix[(x, x)] for x in query_sequence ])

    print("Done")
    print("Etalone score is %d" % etalone_score)
    print("Got %d results, here are top-30 among them:" % len(results))
    print("Score  | Match   | Record")

    for sequence in sorted(results, key=lambda x: x.score, reverse=True)[:30]:
        match = (sequence.score / etalone_score) * 100.0
        print("%6d | %5.3f%% | %s" % (sequence.score, match, sequence.title))

    timer = get_performance_timer()
    for time in [timer.dotplot, timer.regions, timer.align]:
        print(time / cpu_count)

def test_simple():
    #str1 = "TACCGA"
    #str2 = "ACTGAC"

    str1 = "HEAGAWGHEE"
    str2 = "PAWHEAE"

    str1 = "CCATCGCCATCG"
    str2 = "GCATCGGC"

    aligned1, aligned2, score = align(db_seq=str1, query_seq=str2)
    alignment = outline_alignment_for(aligned1, aligned2)
    print(aligned1)
    print(alignment)
    print(aligned2)
    print(score)

def test_rattus():
    library = Fasta(library_path)
    queries = Fasta(query_path)

    query_sequence = str(queries["Rattus"])
    #query_sequence = str(library["sp|P06757|ADH1_RAT Alcohol dehydrogenase 1 OS=Rattus norvegicus GN=Adh1 PE=1 SV=3"])
    #query_sequence = str(library["tr|G0Z9N0|G0Z9N0_MALDO Alcohol dehydrogenase (Fragment) OS=Malus domestica GN=Adh1-1 PE=4 SV=1"])
    #query_sequence = str(library["tr|A0A0J6Z9N4|A0A0J6Z9N4_9MYCO Alcohol dehydrogenase OS=Mycobacterium obuense GN=adh_1 PE=4 SV=1"])
    #query_sequence = str(library["tr|A0A0A9PM57|A0A0A9PM57_ARUDO Adh1 OS=Arundo donax PE=4 SV=1"])
    query_sequence = str(library["tr|Q94462|Q94462_DROBU Alcohol dehydrogenase OS=Drosophila buzzatii GN=Adh1 PE=3 SV=1"])
    #query_sequence = str(library["sp|P48586|ADH1_DROMN Alcohol dehydrogenase 1 OS=Drosophila montana GN=Adh1 PE=3 SV=2"])
    #query_sequence = str(library["tr|W6S386|W6S386_9RHIZ Alcohol dehydrogenase zinc-binding domain protein OS=Rhizobium sp. LPU83 GN=adh1 PE=3 SV=1"])
    #query_sequence = str(library["tr|A0A0E9DQE5|A0A0E9DQE5_CHLTH NADP-dependent alcohol dehydrogenase OS=Chlamydia trachomatis GN=adh_1 PE=3 SV=1"])
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
    process_query()
    #test_simple()
    #test_rattus()

if __name__ == "__main__":
    main()

