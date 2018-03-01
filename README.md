# HMM-Basecaller
MinION, by Oxford Nanopore Technologies, is a third generation sequencer.
Due to its small size and low cost it has become very popular. 
One of its disadvantages in comparison to more expensive and bigger sequencers is a higher error rate in basecalling. 
In order to reduce the error new basecallers are being developed.
Signals received from the MinION are segmented into a list of events and are then translated into a nucleotide sequence by a basecaller. 
One of the basecallers is Nanocall, first open source basecaller. 
Based on the implementation of Nanocall, basecaller described in this work has been made. 
Implemented basecaller relies on Hidden Markov Model to represent the events as the observations and the nucleotides as the hidden states. 
It uses Baum-Welch algorithm to update the transition probabilities 
and Expectation Maximization to update the emission probabilities. 
Having a trained model, Viterbi algorithm is used to decode the sequence of events. 

Use instructions: 
1. git clone https://github.com/anamarijaselak/HMM-Basecaller.git
2. ParseFast5.py should be run in order to preprocess events. Input argument is path to a folder with FAST5 files.
    ParseFast5.py creates new processes files and stores them to folder "ProcessedFiles".
3. Run hmm.cpp
    Input arguments for hmm.cpp are :
        1. path to a directory containing processed files
        2. path to pore model file
        3. path to file Neighbours.txt
        4. path to file BackwardNeighbours.txt
        5. path to file Mask.txt
        6. path to file BackwardMask.txt
        7. path to file indexesForBW
        8. path to a folder where the FASTA files should be stored
