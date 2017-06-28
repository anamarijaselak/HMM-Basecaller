import numpy as np
import h5py
from Bio import pairwise2
import time
from os import listdir
from os.path import isfile, join
from tqdm import tqdm
import sys

'''
N: number of hidden states
'''

def get6mers():
    with open("r73.c.p1.006.ont.model", "r") as ins:
        array = []
        for line in ins:
            array.append(line.split()[0])
    return array

def trainSkipAndStay():
    return 0.1, 0.3

def getTransitionProbability(sixmer_x, sixmer_y, p_stay, p_skip):
    prob = 0

    # pomak za 0 bazi
    if sixmer_x == sixmer_y:
        prob += p_stay
        #return prob

    #pomak za tocno jednu bazu
    if sixmer_x[1:6] == sixmer_y[0:5]:
        prob += (1/4)*(1-p_stay-p_skip)
        #return prob

    #pomak za dvije baze
    if sixmer_x[2:6] == sixmer_y[0:4]:
        prob += (1/16)*(p_skip/(1+p_skip))
        #return prob

    return prob

def defineTransitions():

    transProb = np.zeros((4096, 4096))

    #lista 6-mera, AAAAAA je na indeksu 0
    sixmers = get6mers()

    p_stay, p_skip = trainSkipAndStay()

    for sixmer_x, i in zip(sixmers, range(4096)):
        for sixmer_y, j in zip(sixmers, range(4096)):
            #vjer prijelaza iz sixmer_x u sixmer_y
            transProb[i][j] = getTransitionProbability(sixmer_x, sixmer_y, p_stay, p_skip)

    return transProb

def getNeighbours():
    sixmers = get6mers()

    #iz koga se sve moglo doć u mene
    neighbours = np.zeros((86016), dtype=np.int)
    indexes = np.zeros((86016), dtype=np.int)
    mask = np.ones((86016), dtype=np.int)

    #u koga se sve moze doc iz mene
    backwardNeighbours = np.zeros((86016), dtype=np.int)
    backwardMask = np.ones((86016), dtype=np.int)

    for sixmer, i in zip(sixmers,range(0, 4096)):

        for j in range(0, 21):
            indexes[i*21+j] = i

        #himself
        neighbours[i*21] = i

        #distance 1
        neighbours[i * 21 + 1] = sixmers.index('A' + sixmer[0:5])
        neighbours[i * 21 + 2] = sixmers.index('C' + sixmer[0:5])
        neighbours[i * 21 + 3] = sixmers.index('G'+ sixmer[0:5])
        neighbours[i * 21 + 4] = sixmers.index('T' + sixmer[0:5])

        #distance 2
        neighbours[i * 21 + 5] = sixmers.index('AA'+sixmer[0:4])
        neighbours[i * 21 + 6] = sixmers.index('AC'+sixmer[0:4])
        neighbours[i * 21 + 7] = sixmers.index('AG'+sixmer[0:4])
        neighbours[i * 21 + 8] = sixmers.index('AT'+sixmer[0:4])
        neighbours[i * 21 + 9] = sixmers.index('CA'+sixmer[0:4])
        neighbours[i * 21 + 10] = sixmers.index('CC'+sixmer[0:4])
        neighbours[i * 21 + 11] = sixmers.index('CG'+sixmer[0:4])
        neighbours[i * 21 + 12] = sixmers.index('CT'+sixmer[0:4])
        neighbours[i * 21 + 13] = sixmers.index('GA'+sixmer[0:4])
        neighbours[i * 21 + 14] = sixmers.index('GC'+sixmer[0:4])
        neighbours[i * 21 + 15] = sixmers.index('GG'+sixmer[0:4])
        neighbours[i * 21 + 16] = sixmers.index('GT'+sixmer[0:4])
        neighbours[i * 21 + 17] = sixmers.index('TA'+sixmer[0:4])
        neighbours[i * 21 + 18] = sixmers.index('TC'+sixmer[0:4])
        neighbours[i * 21 + 19] = sixmers.index('TG'+sixmer[0:4])
        neighbours[i * 21 + 20] = sixmers.index('TT'+sixmer[0:4])

    for i in range(0, 4096):
        seen = np.ones((21), dtype=np.int) * (-1)
        for j in range(0, 21):
            if(seen.__contains__(neighbours[i*21 + j])):
                mask[i*21 + j]= 0
            seen[j] = neighbours[i*21 + j]


    #u koga se sve može doć iz sixmer-a
    for sixmer, i in zip(sixmers,range(0, 4096)):

        #himself
        backwardNeighbours[i*21] = i

        #distance 1
        backwardNeighbours[i * 21 + 1] = sixmers.index(sixmer[1:6]+'A')
        backwardNeighbours[i * 21 + 2] = sixmers.index(sixmer[1:6]+'C')
        backwardNeighbours[i * 21 + 3] = sixmers.index(sixmer[1:6]+'G')
        backwardNeighbours[i * 21 + 4] = sixmers.index(sixmer[1:6]+'T')

        #distance 2
        backwardNeighbours[i * 21 + 5] = sixmers.index(sixmer[2:6]+'AA')
        backwardNeighbours[i * 21 + 6] = sixmers.index(sixmer[2:6]+'AC')
        backwardNeighbours[i * 21 + 7] = sixmers.index(sixmer[2:6]+'AG')
        backwardNeighbours[i * 21 + 8] = sixmers.index(sixmer[2:6]+'AT')
        backwardNeighbours[i * 21 + 9] = sixmers.index(sixmer[2:6]+'CA')
        backwardNeighbours[i * 21 + 10] = sixmers.index(sixmer[2:6]+'CC')
        backwardNeighbours[i * 21 + 11] = sixmers.index(sixmer[2:6]+'CG')
        backwardNeighbours[i * 21 + 12] = sixmers.index(sixmer[2:6]+'CT')
        backwardNeighbours[i * 21 + 13] = sixmers.index(sixmer[2:6]+'GA')
        backwardNeighbours[i * 21 + 14] = sixmers.index(sixmer[2:6]+'GC')
        backwardNeighbours[i * 21 + 15] = sixmers.index(sixmer[2:6]+'GG')
        backwardNeighbours[i * 21 + 16] = sixmers.index(sixmer[2:6]+'GT')
        backwardNeighbours[i * 21 + 17] = sixmers.index(sixmer[2:6]+'TA')
        backwardNeighbours[i * 21 + 18] = sixmers.index(sixmer[2:6]+'TC')
        backwardNeighbours[i * 21 + 19] = sixmers.index(sixmer[2:6]+'TG')
        backwardNeighbours[i * 21 + 20] = sixmers.index(sixmer[2:6]+'TT')

    for i in range(0, 4096):
        seen = np.ones((21), dtype=np.int) * (-1)
        for j in range(0, 21):
            if (seen.__contains__(backwardNeighbours[i * 21 + j])):
                backwardMask[i * 21 + j] = 0
            seen[j] = backwardNeighbours[i * 21 + j]

    return neighbours, indexes, mask, backwardNeighbours, backwardMask


def scale_pore_models(models, events):
    #Method of moments applied:
    events_mean = np.sum(events[:, 1]) / len(events[:,1])
    for i in range(0, 4096):
        models[i][0] = events_mean
    return models


def parsePoreModel(filename, events):

    models = np.zeros((4096, 5))
    with open(filename, "r") as ins:
        index = 0
        for line in ins:
            parts = line.split()
            models[index] = parts[1:]
            index += 1

    models = scale_pore_models(models, events)
    return models


def Obs(observation, emission_model,N):
    # racunam vjerojatnost za svaki 6-mer da je emitirao dobiveni observation i vracam to kao numpy array

    sE = time.time()
    mean = observation[1]
    variance = observation[2]

    res = np.zeros(N)

    mean_gaussian = emission_model[:, 0]
    st_deviation_gaussian = emission_model[:, 1]

    mean_inverse_gaussian = emission_model[:, 2]
    st_deviation_inverse_gaussian = emission_model[:, 3]
    lamda_inverse_gaussian = emission_model[:, 4]
    pdf_gaussian = (1/(np.sqrt(2*np.pi*st_deviation_gaussian*st_deviation_gaussian)))*(np.exp((-np.power((mean-mean_gaussian), 2))/(2*st_deviation_gaussian*st_deviation_gaussian)))


    start_time = time.time()
    pdf_inverse_gaussian = np.sqrt((lamda_inverse_gaussian/ (2 * np.pi * np.power(variance, 3)))) * np.exp((
                                                                                                                       -1 * lamda_inverse_gaussian * (
                                                                                                                       variance - mean_inverse_gaussian) * (
                                                                                                                       variance - mean_inverse_gaussian)) / (
                                                                                                                       2 * mean_inverse_gaussian * mean_inverse_gaussian * variance))



    res= pdf_gaussian * pdf_inverse_gaussian
    #print("--- %s Emission ---" % (time.time() - sE))

    return res


class Decoder(object):
    def __init__(self, initialProb, transProb, emissionModel, add):
        self.N = initialProb.shape[1]
        self.initialProb = initialProb
        self.transProb = transProb
        self.emissionModel = emissionModel
        self.neighbours, self.indexes, self.mask, self.backwardN, self.backwardMask = getNeighbours()
        self.add = add
        assert self.initialProb.shape == (1, self.N)
        assert self.transProb.shape == (self.N, self.N)

    def getForwardValues(self, obs):
        print("Start forward step...")
        forwardMatrix = np.ones((self.N, len(obs)))
        forwardMatrix[:, 0] = self.initialProb * Obs(obs[0], self.emissionModel, self.N)

        for o, t in zip(obs[1:], range(1, len(obs))):
            help = forwardMatrix[:, t - 1][self.neighbours] * self.transProb[self.neighbours, self.indexes]*self.mask
            help = help.reshape(4096, 21)
            forwardMatrix[:, t] = np.sum(help, axis=1) * Obs(o, self.emissionModel, self.N)
            forwardMatrix[:, t] = forwardMatrix[:, t] / np.sum(forwardMatrix[:, t])

        #np.savetxt('forwardM2.txt', forwardMatrix, delimiter=',')
        return forwardMatrix


    def getBackwardValues(self, obs):
        print("Start backward step...")
        backwardMatrix = np.ones((self.N, len(obs)))

        for t in range(len(obs) - 1, 0, -1):
            #help = backwardMatrix[:, t][:, np.newaxis] * np.transpose(self.transProb) * Obs(obs[t], self.emissionModel, self.N)[:, np.newaxis]
            #backwardMatrix[:, t - 1] = np.sum(help, axis=0)
            #np.savetxt('backwardM.txt', backwardMatrix, delimiter=',')
            #backwardMatrix[:, t - 1] = backwardMatrix[:, t - 1] / np.sum(backwardMatrix[:, t - 1])

            help = backwardMatrix[:, t][self.backwardN] * self.transProb[self.indexes, self.backwardN] * Obs(obs[t], self.emissionModel, self.N)[self.backwardN]*self.backwardMask
            help = help.reshape(4096, 21)
            backwardMatrix[:, t - 1] = np.sum(help, axis=1)
            backwardMatrix[:, t - 1] = backwardMatrix[:, t - 1] / np.sum(backwardMatrix[:, t - 1])



        return backwardMatrix

    def getPsi(self, t, obs, alfa, beta):
        psi = np.ones((self.N, self.N))
        psi = alfa[:, t][:, np.newaxis] * beta[:, t + 1] * self.transProb * Obs(obs[t + 1], self.emissionModel, self.N)
        sumPsi = np.sum(psi)
        psi = psi / sumPsi

        return psi

    def getGama(self, obs, alfa, beta):
        print("Calculating gama...")
        gama = np.ones((self.N, len(obs) - 1))
        all_psi_sum = np.zeros((self.N, self.N))
        for t in range(0, len(obs) - 1):
            psi_t = self.getPsi(t, obs, alfa, beta)
            all_psi_sum = all_psi_sum + psi_t
            gama[:, t] = np.sum(psi_t, axis=1)

        return gama, all_psi_sum

    def updateParameters(self, obs):
        print("Start updating..")
        gama, psi_summed = self.getGama(obs, self.getForwardValues(obs), self.getBackwardValues(obs))
        gama_summed = np.zeros((self.N, 1))
        gama_summed[:, 0] = np.sum(gama, axis=1)
        self.transProb = psi_summed / gama_summed
        self.initialProb = gama[:, 0]


    def Decode(self, obs):

        p = time.time()

        trellis = np.ones((self.N, len(obs)))
        backpt = np.ones((self.N, len(obs)), 'int32') * -1

        s = time.time()

        help = self.transProb[self.neighbours , self.indexes ] * Obs(obs[0], self.emissionModel, self.N)[self.neighbours] * (1.0/4096.0)
        help = help.reshape(4096, 21)
        trellis[:, 1] = np.amax(help, axis=1)
        axis = np.argmax(help, axis = 1)+self.add
        backpt[:, 1] = self.neighbours[axis]
        trellis[:, 1] = trellis[:, 1] / np.sum(trellis[:, 1])

        for index in range(2, len(obs)):
            help = trellis[:, index - 1][self.neighbours] * self.transProb[self.neighbours, self.indexes] * Obs(obs[index - 1], self.emissionModel, self.N)[self.neighbours]
            help = help.reshape(4096, 21)
            trellis[:, index] = np.amax(help, axis=1)
            axis = np.argmax(help, axis=1) + add
            backpt[:, index] = self.neighbours[axis]
            trellis[:, index] = trellis[:, index] / np.sum(trellis[:, index])

        tokens = []
        #np.savetxt('trellis.txt', trellis[:, 16000], delimiter=',')
        last_state = np.argmax(trellis[:, len(obs) - 1]*Obs(obs[-1], self.emissionModel, self.N))
        tokens.append(last_state)

        for i in reversed(range(1, len(obs))):
            tokens.append(backpt[last_state][i])
            last_state = backpt[last_state][i]

        tokens.reverse()

        t = (time.time() - s) / len(obs)
        #print("--- %s one event ---" % ((t)))

        return tokens

def runViterbi(events, transProb, emission_model, initialProb, add):

    print("Viterbi algorithm started...")
    start_time = time.time()

    d = Decoder(initialProb, transProb, emission_model, add)
    #d.getBackwardValues(events)
    states = []
    nmbOfEvents = len(events)
    chunks = int(nmbOfEvents / 1000) + 1
    for i in tqdm(range(0, chunks)):
        start = i*1000
        end = start + 1000
        if end > nmbOfEvents:
            end = nmbOfEvents
        partial_events = events[start:end:]
        states.extend(d.Decode(partial_events))

    print("Viterbi finished..Output sequence:")
    print("--- %s seconds ---" % (time.time() - start_time))
    return states

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

    print(len(x) / len(y))

    return x, y

def read_fast5(filename):
    'This assumes we are in the right dir.'
    h5 = h5py.File(filename, 'r')

    sampling_rate = h5['UniqueGlobalKey/channel_id'].attrs['sampling_rate']

    reads = h5['Analyses/EventDetection_000/Reads']
    events = np.array(reads[list(reads.keys())[0] + '/Events'])

    basecalled_events = h5['/Analyses/Basecall_1D_000/BaseCalled_template/Events']


    basecalled = np.array(basecalled_events.value[['mean', 'stdv', 'model_state', 'move', 'start',
        'length']])

    # x[i] is feat values for event #i
    x = np.zeros((events.shape[0], 4))

    bcall_idx = 0
    i = 0

    #y[2*i] and y[2*i + 1] are bases for event #i

    y = ['N'] * events.shape[0] * 2
    last_event = 0
    first_event = 0

    for idx, e in enumerate(events):
        if bcall_idx < basecalled.shape[0]:
            b = basecalled[bcall_idx]
            if b[0] == e[2] and b[1] == e[3]: # mean == mean and stdv == stdv
                if bcall_idx == 0:
                    y[2*i-5:2*i] = map(chr,list(b[2])) # initial model state
                    first_event = idx
                bcall_idx += 1
                if b[3] == 1:
                    y[2*i] = chr(b[2][-1])
                if b[3] == 2:
                    y[2*i] = chr(b[2][-2])
                    y[2*i+1] = chr(b[2][-1])
                last_event = idx
        i += 1

    x[:,0] = events['length']
    x[:,1] = events['mean']
    x[:,2] = events['stdv']
    x[:,3] = events['start']

    print(len(x), len(y))

    h5.close()
    return x, y

def getBases(states):

    sixmers = get6mers()
    bases = []
    bases.append(sixmers[states[0]])

    for i in range(1, len(states)):
        before = sixmers[states[i-1]]
        current = sixmers[states[i]]
        if before != current:
            if before[1:6] == current[0:5]:
                bases.append(current[5])
            elif before[2:6] == current[0:4]:
                bases.append(current[4:])
    return "".join(bases)

name, path, poreModelFile = sys.argv

path = "/home/mare/Radna površina/lambda_R9_6_7_16"
poreModelFile = "/home/mare/Radna površina/Projekt/template_median68pA.model"
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

print("Defining transition probabilities...")
transProb = defineTransitions()

initialProb = np.ones((1, 4096))*(1.0/4096.0)

add = np.zeros((4096), dtype=np.int)
for i in range (0, 4096):
    add[i] = i*21

i=0
for file in onlyfiles:
    filename = path + "/"+file
    events, bases = read_fast52(filename)
    emission_model = parsePoreModel(poreModelFile, events)
    states = runViterbi(events, transProb, emission_model, initialProb, add)
    computedBases = getBases(states)
    name = "/data/fastas/Output" + str(i)+".fa"
    text_file = open(name, "w+")
    text_file.write(">%s\n%s\n" % (file, computedBases))
    text_file.close()
    #alignments = pairwise2.align.globalxx(computedBases, real_bases, score_only=True)
    #print(alignments/len(bases))


#filename = "26075_ch100_read112_strand1.fast5"
#events, bases = read_fast5(filename)
#partial_events = events[100:102:]
#partial_bases = bases[200:204]


#states = runViterbi(partial_events, partial_bases)
#bases = getBases(states)
#real_bases = "".join(partial_bases)

#print("Viterbi finished")
#print(bases)
#print("Real bases:")
#print("".join(partial_bases))

#print("Score:")
#alignments = pairwise2.align.globalxx(bases, real_bases, score_only=True)
#print(alignments/len(partial_bases))
