import sys
import numpy as np
import h5py
from os import listdir
from os.path import isfile, join

EPS = 1e-7

def cmp_fp(n1, n2):
    return abs(n1 - n2) < EPS

def read_fast52(filename):
    h5 = h5py.File(filename, 'r')

    sampling_rate = h5['UniqueGlobalKey/channel_id'].attrs['sampling_rate']

    reads = h5['Analyses/EventDetection_000/Reads']
    events = np.array(reads[list(reads.keys())[0] + '/Events'])

    basecalled_events = h5['/Analyses/Basecall_1D_000/BaseCalled_template/Events']

    basecalled = np.array(basecalled_events.value[['mean', 'stdv', 'model_state', 'move', 'start',
                                                   'length']])

    first_basecalled_event = basecalled[0]
    last_basecalled_event = basecalled[-1]

    signal_start = first_basecalled_event['start'] * sampling_rate
    signal_end = last_basecalled_event['start'] * sampling_rate

    start_id = end_id = 0

    for i, e in enumerate(events):
        if cmp_fp(e['start'], signal_start):
            start_id = i
        if cmp_fp(e['start'], signal_end):
            end_id = i
            break

    assert start_id <= end_id

    events = events[start_id:end_id + 1]
    x = np.zeros((events.shape[0], 4))
    x[:, 0] = events['length']
    x[:, 1] = events['mean']
    x[:, 2] = events['stdv']
    x[:, 3] = events['start']

    bcall_idx = 0
    i = 0

    y = ['N'] * events.shape[0] * 2

    for e in events:
        if bcall_idx < basecalled.shape[0]:
            b = basecalled[bcall_idx]
            if b[0] == e[2] and b[1] == e[3]:  # mean == mean and stdv == stdv
                assert cmp_fp(b['start'] * sampling_rate, e['start'])
                if bcall_idx == 0:
                    y[2 * i - 5:2 * i] = map(chr, list(b[2]))  # initial model state
                bcall_idx += 1
                if b[3] == 1:
                    y[2 * i] = chr(b[2][-1])
                if b[3] == 2:
                    y[2 * i] = chr(b[2][-2])
                    y[2 * i + 1] = chr(b[2][-1])
        i += 1

    h5.close()

    while "N" in y: y.remove("N")

   # print(len(x) / len(y))

    return x, y

path = sys.argv[1]
print(path)
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
i=0
for file in onlyfiles:
    #1 file = 1 read
    i = i+1
    print("File")
    print(i)
    print("Start processing the file...")
    filename = path + "/"+file
    events, bases = read_fast52(filename)
    f = open(path + "/ProcessedFiles/" +file, "w")
    f.write(str(events.size / 4))
    f.write("\n")
    np.savetxt(path + "/ProcessedFiles/" +file, events)