## Encrypt and decrypt strings to DNA Sequences

### SYNOPSIS

```perl
use Crypt::DNASequence;
Crypt::DNASequence->encrypt($input_filename);
Crypt::DNASequence->decrypt($encrypted_file);
```

### DESCRIPTION

The module is naive and just for fun. It transforms text strings into DNA sequences. A DNA sequence
is composed of four nucleotides which are represented as A, T, C, G. If we transform 
"abcdefghijklmnopqistuvwxyzABCDEFGHIJKLMNOPQISTUVWXYZ", the corresponding sequence would be:

    TAGACCTGTCGGTGGGTGGTTCCCCGTATGAAGGGGGGGGCTATGAGTCTACAGACTGACAGGGGGAGTC
    AGCGAGTTAGCTAGTAAGCAAGTCCGCCCGTGCGCGCGTTCGCTCGTACGCACGTCGTCCGTTGCGCGGT
    CGGTCTGTTAGTCAGTTCTTCCTTTGTTCGACGTACAGACGTACATACGAACAAACGCCCACCCGGCCAG
    CCGTCCATCCGACCAACCGCGGCCGGTGCCAGGGCGGGCTGGTAGGCAGGTCTGCCTGTGTGCGGGGTGG
    GGCCCACGCCGCACAGCAGCAGGGGGGGGGGGGGGC

or

    CGAGTTCACTAACAAACAACCTTTTACGCAGGAAAAAAAATCGCAGACTCGTGAGTCAGTGAAAAAGACT
    GATAGACCGATCGACGGATGGACTTATTTACATATATACCTATCTACGTATGTACTACTTACCATATAAC
    TAACTCACCGACTGACCTCCTTCCCACCTAGTACGTGAGTACGTGCGTAGGTGGGTATTTGTTTAATTGA
    TTACTTGCTTAGTTGGTTATAATTAACATTGAAATAAATCAACGAATGAACTCATTCACACATAAAACAA
    AATTTGTATTATGTGATGATGAAAAAAAAAAAAAAT
  
The transformation is not unique due to a random mapping, but all the transformed sequences can be 
decrypted correctly to the origin string.

### ALGORITHM

First, text file are compressed into gzip files when encrypting and DNA sequences
are first decrypted into gzip files and then uncompressed into normal text file.

The algorithm behind the module is simple. Two binary bits are used to represent a nucleotide such as '00' for A, '01' for C. 
If you have some knowledge of molecular biology, you would know that A only matches to T and C only matches to G.
So if '00' is chosen to be A, then '11' should be used to represent 'T'. In the module, the correspondence between binary bits
and nucleotides are applied randomly. The information of the correspondence dictionary is also stored in the finnal sequence.

Here is the procedure for encryption. 1. Split a string into a set of letters or charactors. 2. For each letter, convert to
its binary form and transform to ATCG every two bits using a randomly generated dictionary. The dictionary may looks like:

```perl
$dict = { '00' => 'A',
          '11' => 'T',
          '01' => 'C',
          '10' => 'G' };
```

3. Join the A, T, G, C as a single sequence. 4. Find the first nucleotide of the sequence. 5, Find the number of the first nucleotide
in the sequence. 6. There is a database storing all arrangements of '00', '11', '01', '10'. 7. Calculate the index value from
the number of the first nucleotide by mod calculation. 8. Retrieve the arrangement with the index value, map them to the dictionary and get four nucleotides. E.g. the first nucleotide of the sequence is G. The number of G in the sequence is 40. The number of all arrangement
in the database is 24. Then we calculate the index value by 40 % 24 = 16. Then the 16th arrangement is retrieved and may looks like
['01', '11', '10', '00']. The four items in the array are mapped to the dictionary to be four nucleotides such as CTGA. Note this information
can be used in the decryption procedure.
9. Put the first two nucleotides at the begining of the sequence and the last two nucleotides at the end of the sequence. 10. That
is the finnal seuqence.

Here is the procedure for decryption. 1. Extract the first two and the last two nucleotides fromt the sequence. E.g. CT and GA. 
2. Count the number of the first nucleotide in the real sequence, e.g., 40 for G. 3. Use this number to calculate the index in the 
arrangement database, e.g., 16. 4. find the dictionary, i.e. a dictionary is generated from the 16th arrangement ['01', '11', '10', '00'] and
CTGA. 5. Translate the DNA sequence according the dictionary into binary bit form and finnaly to the orgin format.

### Subroutines

`Crypt::DNASequence->encrypt($input_file)`: encrypt the text file to DNA sequence

`Crypt::DNASequence->decrypt($encrypted_file)`: decrypt the DNA sequence to the origin text file.
        