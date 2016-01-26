# Introduction

overall 76 parameters generate from previous research work.

# Details

An example of each parameter is attached in the next section.

|Name|Description|
|:---|:-------------------------------------------------------------------------------------------| 
|MIID|The ID of miRNA, usually like this "three letters of species abbreciation-MIR-serial number"|
|PRIID|The ID of pri-miRNA, usually like this "three letters of species abbreciation-mir-serial number"| 
|MISEQ|The sequence of miRNA | 
|MISTR_1|The first strand of miRNA duplex secondary structure | 
|MISTR_2|The second strand of miRNA duplex secondary structure | |MISTR_3|The third strand of miRNA duplex secondary structure | 
|MISTR_4|The fourth strand of miRNA duplex secondary structure | 
|PRESEQ|The sequence of miRNA precusors, started with miRNA duplex | |PRESTR_1|The first strand of pre-miRNA secondary structure | 
|PRESTR_2|The second strand of pre-miRNA secondary structure | 
|PRESTR_3|The third strand of pre-miRNA secondary structure | 
|PRESTR_4|The fourth strand of pre-miRNA secondary structure | 
|PRISEQ|The sequence of primary miRNA | 
|PRISTR_1|The first strand of pri-miRNA secondary structure | 
|PRISTR_2|The second strand of pri-miRNA secondary structure | 
|PRISTR_3|The third strand of pri-miRNA secondary structure | 
|PRISTR_4|The fourth strand of pri-miRNA secondary structure | 
|MILENGTH|The length of miRNA (nt) | 
|PRELENGTH|The length of pre-miRNA (nt) | 
|PRILENGTH|The length of pri-miRNA (nt) | 
|LENGTH_BASALSEGMENT|The size of The Basal Segment secondary structure (bp), see our paper for more details. | 
|LENGTH_LOWERSTEM|The size of The Lower Stem secondary structure (bp), see our paper for more details. | 
|LENGTH_UPPERSTEM|The size of The Upper Stem secondary structure (bp), see our paper for more details. | 
|LENGTH_TOPSTEM|The size of The Top Stem secondary structure (bp), see our paper for more details. | 
|LENGTH_TERMINALLOOP|The size of The Terminal Loop (nt), see our paper for more details. | 
|MIPAIRS|The number of nucleotide pairs in miRNA duplex (bp) | 
|PREPAIRS|The number of nucleotide pairs in pre-miRNA secondary structure (bp) | 
|PRIPAIRS|The number of nucleotide pairs in pri-miRNA secondary structure (bp) | 
|PREMFE|The Minimal Free Energy of pre-miRNA (kcal/mol) | 
|PRIMFE|The Minimal Free Energy of pri-miRNA (kcal/mol) | 
|MIGC|The GC Content of miRNA | |PREGC|The GC Content of pre-miRNA | 
|PRIGC|The GC Content of pri-miRNA | 
|MINTCONTENT_A|The Nucleotide Content of "A" in miRNA | 
|MINTCONTENT_C|The Nucleotide Content of "C" in miRNA | 
|MINTCONTENT_G|The Nucleotide Content of "G" in miRNA | 
|MINTCONTENT_U|The Nucleotide Content of "U" in miRNA | 
|PRENTCONTENT_A|The Nucleotide Content of "A" in pre-miRNA | 
|PRENTCONTENT_C|The Nucleotide Content of "C" in pre-miRNA | 
|PRENTCONTENT_G|The Nucleotide Content of "G" in pre-miRNA | 
|PRENTCONTENT_U|The Nucleotide Content of "U" in pre-miRNA | 
|PRINTCONTENT_A|The Nucleotide Content of "A" in pri-miRNA | 
|PRINTCONTENT_C|The Nucleotide Content of "C" in pri-miRNA | 
|PRINTCONTENT_G|The Nucleotide Content of "G" in pri-miRNA | 
|PRINTCONTENT_U|The Nucleotide Content of "U" in pri-miRNA | 
|MIINTERNALLOOP|The size of the biggest internal loop in miRNA duplex (nt) | 
|PREINTERNALLOOP|The size of the biggest internal loop in pre-miRNA secondary structure (nt) | 
|PRIINTERNALLOOP|The size of the biggest internal loop in pri-miRNA secondary structure (nt) | 
|INTERNALLOOP_LOWERSTEM|The size of the biggest internal loop in the Lower Stem secondary structure (nt) | 
|INTERNALLOOP_TOPSTEM|The size of the biggest internal loop in the Top Stem secondary structure (nt) | 
|MIINTERNALLOOPNUMBER|The number of internal loop in miRNA duplex | 
|PREINTERNALLOOPNUMBER|The number of internal loop in pre-miRNA secondary structure | 
|PRIINTERNALLOOPNUMBER|The number of internal loop in pri-miRNA secondary structure | 
|INTERNALLOOPNUMBER_LOWERSTEM|The number of internal loop in the Lower Stem secondary structure | 
|INTERNALLOOPNUMBER_TOPSTEM|The number of internal loop in the Top Stem secondary structure | 
|MIUNPAIREDBASES|The total number of unpaired bases in miRNA duplex (nt) | 
|PREUNPAIREDBASES|The total number of unpaired bases in pre-miRNA secondary structure (nt) | 
|PRIUNPAIREDBASES|The total number of unpaired bases in pri-miRNA secondary structure (nt) | 
|UNPAIREDBASES_LOWERSTEM|The total number of unpaired bases in the Lower Stem secondary structure (nt) | 
|UNPAIREDBASES_TOPSTEM|The total number of unpaired bases in the Top Stem secondary structure (nt) | 
|MIUNPAIREDRATE|The unpaired rate in miRNA duplex | 
|PREUNPAIREDRATE|The unpaired rate in pre-miRNA secondary structure | 
|PRIUNPAIREDRATE|The unpaired rate in pri-miRNA secondary structure | 
|UNPAIREDRATE_LOWERSTEM|The unpaired rate in the Lower Stem secondary structure | 
|UNPAIREDRATE_TOPSTEM|The unpaired rate in the Top Stem secondary structure | 
|MIGU|The number of G:U wobbles in miRNA duplex | 
|PREGU|The number of G:U wobbles in pre-miRNA secondary structure | 
|PRIGU|The number of G:U wobbles in pri-miRNA secondary structure | 
|STRAND|The strand that miRNA located in pri-miRNA, it could either be 5' or 3' | |FIRSTBASE|The first nucleotide of miRNA | 
|PENULTIMATEPOSITION|The second last nucleotide of 3' overhang | 
|TERMINALNUCLEOTIDE|The last nucleotide of 3' overhang | |MISTART|The start position of miRNA in pri-miRNA sequences | 
|MIEND|The end position of miRNA IN pri-miRNA sequences | 
|UPPERSTART|The start position of miRNA duplex in pri-miRNA secondary strucutre | 
|UPPEREND|The end position of miRNA duplex in pri-miRNA secondary strucutre | 
|STABILITY|The stability of two terminal of miRNA duplex, it calculated from the rate of total number of hydrogen bonds|

# Example

|Name:|Example| 
|:----|:--------------| 
|MIID:|cel-MIR-105_22| 
|PRIID:|cel_1-99-mir-1| 
|MISEQ:|ugagguaguagguuguauaguu| 
|MISTR_1:|U GG | 
|MISTR_2:|UGAGGUAG A UUGUAUAGUU| 
|MISTR_3:|auuccauc u aacguaucaa| 
|MISTR_4:|u uu | 
|PRESEQ:|ugagguaguagguuguauaguuuggaauauuaccaccggugaacuaugcaauuuucuaccuua| 
|PRESTR_1:|U GG --- aaua | 
|PRESTR_2:|UGAGGUAG A UUGUAUAGUUu gg u| 
|PRESTR_3:|auuccauc u aacguaucaag cc u| 
|PRESTR_4:|u uu ugg acca | 
|PRISEQ:|uacacuguggauccggugagguaguagguuguauaguuuggaauauuaccaccggugaacuaugcaauuuucuaccuuaccggagacagaacucuucga| 
|PRISTR_1:|uacacug---- ga U GG --- aaua | 
|PRISTR_2:|ug uccggUGAGGUAG A UUGUAUAGUUu gg u| 
|PRISTR_3:|ac aggccauuccauc u aacguaucaag cc u| 
|PRISTR_4:|agcuucucaag ag u uu ugg acca | 
|MILENGTH:|22 | 
|PRELENGTH:|63 | 
|PRILENGTH:|99 | 
|LENGTH_BASALSEGMENT:|11 | 
|LENGTH_LOWERSTEM:|9 | 
|LENGTH_UPPERSTEM:|22 | 
|LENGTH_TOPSTEM:|6 | 
|LENGTH_TERMINALLOOP:|10 | 
|MIPAIRS:|19 | 
|PREPAIRS:|22 | 
|PRIPAIRS:|29 | 
|PREMFE:|-25.45 | 
|PRIMFE:|-40 | 
|MIGC:|0.3636 | 
|PREGC:|0.3651 | 
|PRIGC:|0.4343 | 
|MINTCONTENT_A:|0.2273 | 
|MINTCONTENT_C:|0.0000 | 
|MINTCONTENT_G:|0.3636 | 
|MINTCONTENT_U:|0.4091 | 
|PRENTCONTENT_A:|0.2698 | 
|PRENTCONTENT_C:|0.1429 | 
|PRENTCONTENT_G:|0.2222 | 
|PRENTCONTENT_U:|0.3651 | 
|PRINTCONTENT_A:|0.2626 | 
|PRINTCONTENT_C:|0.1919 | 
|PRINTCONTENT_G:|0.2424 | 
|PRINTCONTENT_U:|0.3030 | 
|MIINTERNALLOOP:|2 | 
|PREINTERNALLOOP:|3 | 
|PRIINTERNALLOOP:|3 | 
|INTERNALLOOP_LOWERSTEM:|2 | 
|INTERNALLOOP_TOPSTEM:|3 | 
|MIINTERNALLOOPNUMBER:|2 | 
|PREINTERNALLOOPNUMBER:|2 | 
|PRIINTERNALLOOPNUMBER:|4 | 
|INTERNALLOOPNUMBER_LOWERSTEM:|1 | 
|INTERNALLOOPNUMBER_TOPSTEM:|1 | 
|MIUNPAIREDBASES:|6 | 
|PREUNPAIREDBASES:|9 | 
|PRIUNPAIREDBASES:|13 | 
|UNPAIREDBASES_LOWERSTEM:|4 | 
|UNPAIREDBASES_TOPSTEM:|3 | 
|MIUNPAIREDRATE:|0.1364 | 
|PREUNPAIREDRATE:|0.1698 | 
|PRIUNPAIREDRATE:|0.1831 | 
|UNPAIREDRATE_LOWERSTEM:|0.2222 | 
|UNPAIREDRATE_TOPSTEM:|0.3333 | 
|MIGU:|2 | 
|PREGU:|3 | 
|PRIGU:|3 | 
|STRAND:|5 | 
|FIRSTBASE:|U | 
|PENULTIMATEPOSITION:|G | 
|TERMINALNUCLEOTIDE:|U | 
|MISTART:|17 | 
|MIEND:|38 | 
|UPPERSTART:|21 | 
|UPPEREND:|42 | 
|STABILITY:|1.00 |
