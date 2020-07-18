**createTargetingModel** - *Creates a TargetingModel*

Description
--------------------

`createTargetingModel` creates a 5-mer `TargetingModel`.


Usage
--------------------
```
createTargetingModel(
db,
model = c("s", "rs"),
sequenceColumn = "sequence_alignment",
germlineColumn = "germline_alignment_d_mask",
vCallColumn = "v_call",
multipleMutation = c("independent", "ignore"),
minNumMutations = 50,
minNumSeqMutations = 500,
modelName = "",
modelDescription = "",
modelSpecies = "",
modelCitation = "",
modelDate = NULL
)
```

Arguments
-------------------

db
:   data.frame containing sequence data.

model
:   type of model to create. The default model, "s", 
builds a model by counting only silent mutations. `model="s"`
should be used for data that includes functional sequences.
Setting `model="rs"` creates a model by counting both 
replacement and silent mutations and may be used on fully 
non-functional sequence data sets.

sequenceColumn
:   name of the column containing IMGT-gapped sample sequences.

germlineColumn
:   name of the column containing IMGT-gapped germline sequences.

vCallColumn
:   name of the column containing the V-segment allele calls.

multipleMutation
:   string specifying how to handle multiple mutations occuring 
within the same 5-mer. If `"independent"` then multiple 
mutations within the same 5-mer are counted indepedently. 
If `"ignore"` then 5-mers with multiple mutations are 
excluded from the otal mutation tally.

minNumMutations
:   minimum number of mutations required to compute the 5-mer 
substitution rates. If the number of mutations for a 5-mer
is below this threshold, its substitution rates will be 
estimated from neighboring 5-mers. Default is 50.

minNumSeqMutations
:   minimum number of mutations in sequences containing each 5-mer
to compute the mutability rates. If the number is smaller 
than this threshold, the mutability for the 5-mer will be 
inferred. Default is 500.

modelName
:   name of the model.

modelDescription
:   description of the model and its source data.

modelSpecies
:   genus and species of the source sequencing data.

modelCitation
:   publication source.

modelDate
:   date the model was built. If `NULL` the current date
will be used.




Value
-------------------

A [TargetingModel](TargetingModel-class.md) object.


Details
-------------------

**Caution: The targeting model functions do NOT support ambiguous 
characters in their inputs. You MUST make sure that your input and germline
sequences do NOT contain ambiguous characters (especially if they are
clonal consensuses returned from `collapseClones`).**


References
-------------------


1. Yaari G, et al. Models of somatic hypermutation targeting and substitution based
on synonymous mutations from high-throughput immunoglobulin sequencing data.
Front Immunol. 2013 4(November):358.
 



Examples
-------------------

```R
# Subset example data to one isotype and sample as a demo
data(ExampleDb, package="alakazam")
db <- subset(ExampleDb, c_call == "IGHA" & sample_id == "-1h")

# Create model using only silent mutations and ignore multiple mutations
model <- createTargetingModel(db, model="s", sequenceColumn="sequence_alignment",
germlineColumn="germline_alignment_d_mask",
vCallColumn="v_call", multipleMutation="ignore")

```

*Warning*:Insufficient number of mutations to infer some 5-mers. Filled with 0. 
```R


# Access and view mutability estimates (not run)
print(model@mutability)

```


```
       AAAAA        AAAAC        AAAAG        AAAAT        AAAAN        AAACA        AAACC        AAACG        AAACT        AAACN        AAAGA        AAAGC        AAAGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAAGT        AAAGN        AAATA        AAATC        AAATG        AAATT        AAATN        AAANA        AAANC        AAANG        AAANT        AAANN        AACAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACAC        AACAG        AACAT        AACAN        AACCA        AACCC        AACCG        AACCT        AACCN        AACGA        AACGC        AACGG        AACGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACGN        AACTA        AACTC        AACTG        AACTT        AACTN        AACNA        AACNC        AACNG        AACNT        AACNN        AAGAA        AAGAC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0027781351 0.0000000000 
       AAGAG        AAGAT        AAGAN        AAGCA        AAGCC        AAGCG        AAGCT        AAGCN        AAGGA        AAGGC        AAGGG        AAGGT        AAGGN 
0.0023072876 0.0000000000 0.0012713557 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0002851050 0.0002066547 
       AAGTA        AAGTC        AAGTG        AAGTT        AAGTN        AAGNA        AAGNC        AAGNG        AAGNT        AAGNN        AATAA        AATAC        AATAG 
0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0051708243 0.0000000000 0.0005768219 0.0021539494 0.0019753989 0.0000000000 0.0000000000 0.0000000000 
       AATAT        AATAN        AATCA        AATCC        AATCG        AATCT        AATCN        AATGA        AATGC        AATGG        AATGT        AATGN        AATTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005803880 
       AATTC        AATTG        AATTT        AATTN        AATNA        AATNC        AATNG        AATNT        AATNN        AANAA        AANAC        AANAG        AANAT 
0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0001450970           NA           NA           NA           NA 
       AANAN        AANCA        AANCC        AANCG        AANCT        AANCN        AANGA        AANGC        AANGG        AANGT        AANGN        AANTA        AANTC 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       AANTG        AANTT        AANTN        AANNA        AANNC        AANNG        AANNT        AANNN        ACAAA        ACAAC        ACAAG        ACAAT        ACAAN 
          NA           NA           NA           NA           NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACACA        ACACC        ACACG        ACACT        ACACN        ACAGA        ACAGC        ACAGG        ACAGT        ACAGN        ACATA        ACATC        ACATG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0000000000 0.0000000000 0.0000000000 
       ACATT        ACATN        ACANA        ACANC        ACANG        ACANT        ACANN        ACCAA        ACCAC        ACCAG        ACCAT        ACCAN        ACCCA 
0.0000000000 0.0000000000 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0043115900 0.0043115900 0.0043115900 0.0012186597 0.0035383574 0.0027651248 
       ACCCC        ACCCG        ACCCT        ACCCN        ACCGA        ACCGC        ACCGG        ACCGT        ACCGN        ACCTA        ACCTC        ACCTG        ACCTT 
0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0012186597 0.0012186597 0.0012186597 0.0043115900 
       ACCTN        ACCNA        ACCNC        ACCNG        ACCNT        ACCNN        ACGAA        ACGAC        ACGAG        ACGAT        ACGAN        ACGCA        ACGCC 
0.0019918923 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0009606561 0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0025595886 0.0000000000 
       ACGCG        ACGCT        ACGCN        ACGGA        ACGGC        ACGGG        ACGGT        ACGGN        ACGTA        ACGTC        ACGTG        ACGTT        ACGTN 
0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 
       ACGNA        ACGNC        ACGNG        ACGNT        ACGNN        ACTAA        ACTAC        ACTAG        ACTAT        ACTAN        ACTCA        ACTCC        ACTCG 
0.0047164545 0.0000000000 0.0005768219 0.0024290355 0.0019305780 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACTCT        ACTCN        ACTGA        ACTGC        ACTGG        ACTGT        ACTGN        ACTTA        ACTTC        ACTTG        ACTTT        ACTTN        ACTNA 
0.0000000000 0.0000000000 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005389060 
       ACTNC        ACTNG        ACTNT        ACTNN        ACNAA        ACNAC        ACNAG        ACNAT        ACNAN        ACNCA        ACNCC        ACNCG        ACNCT 
0.0005389060 0.0005389060 0.0005389060 0.0005389060           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ACNCN        ACNGA        ACNGC        ACNGG        ACNGT        ACNGN        ACNTA        ACNTC        ACNTG        ACNTT        ACNTN        ACNNA        ACNNC 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ACNNG        ACNNT        ACNNN        AGAAA        AGAAC        AGAAG        AGAAT        AGAAN        AGACA        AGACC        AGACG        AGACT        AGACN 
          NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AGAGA        AGAGC        AGAGG        AGAGT        AGAGN        AGATA        AGATC        AGATG        AGATT        AGATN        AGANA        AGANC        AGANG 
0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0001563454 0.0001563454 0.0001563454 
       AGANT        AGANN        AGCAA        AGCAC        AGCAG        AGCAT        AGCAN        AGCCA        AGCCC        AGCCG        AGCCT        AGCCN        AGCGA 
0.0001563454 0.0001563454 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCGC        AGCGG        AGCGT        AGCGN        AGCTA        AGCTC        AGCTG        AGCTT        AGCTN        AGCNA        AGCNC        AGCNG        AGCNT 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCNN        AGGAA        AGGAC        AGGAG        AGGAT        AGGAN        AGGCA        AGGCC        AGGCG        AGGCT        AGGCN        AGGGA        AGGGC 
0.0073813816 0.0009606561 0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 
       AGGGG        AGGGT        AGGGN        AGGTA        AGGTC        AGGTG        AGGTT        AGGTN        AGGNA        AGGNC        AGGNG        AGGNT        AGGNN 
0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0047164545 0.0000000000 0.0005768219 0.0024290355 0.0019305780 
       AGTAA        AGTAC        AGTAG        AGTAT        AGTAN        AGTCA        AGTCC        AGTCG        AGTCT        AGTCN        AGTGA        AGTGC        AGTGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0016522853 0.0016522853 0.0016522853 
       AGTGT        AGTGN        AGTTA        AGTTC        AGTTG        AGTTT        AGTTN        AGTNA        AGTNC        AGTNG        AGTNT        AGTNN        AGNAA 
0.0016522853 0.0016522853 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0004130713 0.0004130713 0.0004130713 0.0004130713 0.0004130713           NA 
       AGNAC        AGNAG        AGNAT        AGNAN        AGNCA        AGNCC        AGNCG        AGNCT        AGNCN        AGNGA        AGNGC        AGNGG        AGNGT 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       AGNGN        AGNTA        AGNTC        AGNTG        AGNTT        AGNTN        AGNNA        AGNNC        AGNNG        AGNNT        AGNNN        ATAAA        ATAAC 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 0.0000000000 0.0000000000 
       ATAAG        ATAAT        ATAAN        ATACA        ATACC        ATACG        ATACT        ATACN        ATAGA        ATAGC        ATAGG        ATAGT        ATAGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATATA        ATATC        ATATG        ATATT        ATATN        ATANA        ATANC        ATANG        ATANT        ATANN        ATCAA        ATCAC        ATCAG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCAT        ATCAN        ATCCA        ATCCC        ATCCG        ATCCT        ATCCN        ATCGA        ATCGC        ATCGG        ATCGT        ATCGN        ATCTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCTC        ATCTG        ATCTT        ATCTN        ATCNA        ATCNC        ATCNG        ATCNT        ATCNN        ATGAA        ATGAC        ATGAG        ATGAT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 
       ATGAN        ATGCA        ATGCC        ATGCG        ATGCT        ATGCN        ATGGA        ATGGC        ATGGG        ATGGT        ATGGN        ATGTA        ATGTC 
0.0010441708 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0024857935 0.0007568269 0.0148040595 0.0000000000 
       ATGTG        ATGTT        ATGTN        ATGNA        ATGNC        ATGNG        ATGNT        ATGNN        ATTAA        ATTAC        ATTAG        ATTAT        ATTAN 
0.0000000000 0.0000000000 0.0037010149 0.0049436394 0.0000000000 0.0005768219 0.0027041215 0.0020561457 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTCA        ATTCC        ATTCG        ATTCT        ATTCN        ATTGA        ATTGC        ATTGG        ATTGT        ATTGN        ATTTA        ATTTC        ATTTG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTTT        ATTTN        ATTNA        ATTNC        ATTNG        ATTNT        ATTNN        ATNAA        ATNAC        ATNAG        ATNAT        ATNAN        ATNCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000           NA           NA           NA           NA           NA           NA 
       ATNCC        ATNCG        ATNCT        ATNCN        ATNGA        ATNGC        ATNGG        ATNGT        ATNGN        ATNTA        ATNTC        ATNTG        ATNTT 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ATNTN        ATNNA        ATNNC        ATNNG        ATNNT        ATNNN        ANAAA        ANAAC        ANAAG        ANAAT        ANAAN        ANACA        ANACC 
          NA           NA           NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANACG        ANACT        ANACN        ANAGA        ANAGC        ANAGG        ANAGT        ANAGN        ANATA        ANATC        ANATG        ANATT        ANATN 
0.0000000000 0.0000000000 0.0000000000 0.0010554064 0.0010554064 0.0010554064 0.0010554064 0.0010554064 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANANA        ANANC        ANANG        ANANT        ANANN        ANCAA        ANCAC        ANCAG        ANCAT        ANCAN        ANCCA        ANCCC        ANCCG 
0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0029232429 0.0029232429 0.0029232429 0.0021500103 0.0027299347 0.0025366266 0.0025366266 0.0025366266 
       ANCCT        ANCCN        ANCGA        ANCGC        ANCGG        ANCGT        ANCGN        ANCTA        ANCTC        ANCTG        ANCTT        ANCTN        ANCNA 
0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0021500103 0.0021500103 0.0021500103 0.0029232429 0.0023433185 0.0025366266 
       ANCNC        ANCNG        ANCNT        ANCNN        ANGAA        ANGAC        ANGAG        ANGAT        ANGAN        ANGCA        ANGCC        ANGCG        ANGCT 
0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0016422107 0.0000000000 0.0023072876 0.0000000000 0.0009873746 0.0025595886 0.0000000000 0.0000000000 0.0083306927 
       ANGCN        ANGGA        ANGGC        ANGGG        ANGGT        ANGGN        ANGTA        ANGTC        ANGTG        ANGTT        ANGTN        ANGNA        ANGNC 
0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0048868432 0.0000000000 
       ANGNG        ANGNT        ANGNN        ANTAA        ANTAC        ANTAG        ANTAT        ANTAN        ANTCA        ANTCC        ANTCG        ANTCT        ANTCN 
0.0005768219 0.0024290355 0.0019731751 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANTGA        ANTGC        ANTGG        ANTGT        ANTGN        ANTTA        ANTTC        ANTTG        ANTTT        ANTTN        ANTNA        ANTNC        ANTNG 
0.0009519773 0.0009519773 0.0009519773 0.0009519773 0.0009519773 0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0002742686 0.0002742686 0.0002742686 
       ANTNT        ANTNN        ANNAA        ANNAC        ANNAG        ANNAT        ANNAN        ANNCA        ANNCC        ANNCG        ANNCT        ANNCN        ANNGA 
0.0002742686 0.0002742686           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ANNGC        ANNGG        ANNGT        ANNGN        ANNTA        ANNTC        ANNTG        ANNTT        ANNTN        ANNNA        ANNNC        ANNNG        ANNNT 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ANNNN        CAAAA        CAAAC        CAAAG        CAAAT        CAAAN        CAACA        CAACC        CAACG        CAACT        CAACN        CAAGA        CAAGC 
          NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAAGG        CAAGT        CAAGN        CAATA        CAATC        CAATG        CAATT        CAATN        CAANA        CAANC        CAANG        CAANT        CAANN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACAA        CACAC        CACAG        CACAT        CACAN        CACCA        CACCC        CACCG        CACCT        CACCN        CACGA        CACGC        CACGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACGT        CACGN        CACTA        CACTC        CACTG        CACTT        CACTN        CACNA        CACNC        CACNG        CACNT        CACNN        CAGAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0009606561 
       CAGAC        CAGAG        CAGAT        CAGAN        CAGCA        CAGCC        CAGCG        CAGCT        CAGCN        CAGGA        CAGGC        CAGGG        CAGGT 
0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 0.0000000000 0.0000000000 0.0024857935 
       CAGGN        CAGTA        CAGTC        CAGTG        CAGTT        CAGTN        CAGNA        CAGNC        CAGNG        CAGNT        CAGNN        CATAA        CATAC 
0.0007568269 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0046447092 0.0000000000 0.0005768219 0.0027041215 0.0019814132 0.0000000000 0.0000000000 
       CATAG        CATAT        CATAN        CATCA        CATCC        CATCG        CATCT        CATCN        CATGA        CATGC        CATGG        CATGT        CATGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CATTA        CATTC        CATTG        CATTT        CATTN        CATNA        CATNC        CATNG        CATNT        CATNN        CANAA        CANAC        CANAG 
0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0001450970           NA           NA           NA 
       CANAT        CANAN        CANCA        CANCC        CANCG        CANCT        CANCN        CANGA        CANGC        CANGG        CANGT        CANGN        CANTA 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       CANTC        CANTG        CANTT        CANTN        CANNA        CANNC        CANNG        CANNT        CANNN        CCAAA        CCAAC        CCAAG        CCAAT 
          NA           NA           NA           NA           NA           NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCAAN        CCACA        CCACC        CCACG        CCACT        CCACN        CCAGA        CCAGC        CCAGG        CCAGT        CCAGN        CCATA        CCATC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0000000000 0.0000000000 
       CCATG        CCATT        CCATN        CCANA        CCANC        CCANG        CCANT        CCANN        CCCAA        CCCAC        CCCAG        CCCAT        CCCAN 
0.0000000000 0.0000000000 0.0000000000 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCCA        CCCCC        CCCCG        CCCCT        CCCCN        CCCGA        CCCGC        CCCGG        CCCGT        CCCGN        CCCTA        CCCTC        CCCTG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCTT        CCCTN        CCCNA        CCCNC        CCCNG        CCCNT        CCCNN        CCGAA        CCGAC        CCGAG        CCGAT        CCGAN        CCGCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 0.0010441708 0.0022726072 
       CCGCC        CCGCG        CCGCT        CCGCN        CCGGA        CCGGC        CCGGG        CCGGT        CCGGN        CCGTA        CCGTC        CCGTG        CCGTT 
0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 
       CCGTN        CCGNA        CCGNC        CCGNG        CCGNT        CCGNN        CCTAA        CCTAC        CCTAG        CCTAT        CCTAN        CCTCA        CCTCC 
0.0037010149 0.0048718941 0.0000000000 0.0005768219 0.0024290355 0.0019694379 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCTCG        CCTCT        CCTCN        CCTGA        CCTGC        CCTGG        CCTGT        CCTGN        CCTTA        CCTTC        CCTTG        CCTTT        CCTTN 
0.0000000000 0.0000000000 0.0000000000 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCTNA        CCTNC        CCTNG        CCTNT        CCTNN        CCNAA        CCNAC        CCNAG        CCNAT        CCNAN        CCNCA        CCNCC        CCNCG 
0.0005389060 0.0005389060 0.0005389060 0.0005389060 0.0005389060           NA           NA           NA           NA           NA           NA           NA           NA 
       CCNCT        CCNCN        CCNGA        CCNGC        CCNGG        CCNGT        CCNGN        CCNTA        CCNTC        CCNTG        CCNTT        CCNTN        CCNNA 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       CCNNC        CCNNG        CCNNT        CCNNN        CGAAA        CGAAC        CGAAG        CGAAT        CGAAN        CGACA        CGACC        CGACG        CGACT 
          NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGACN        CGAGA        CGAGC        CGAGG        CGAGT        CGAGN        CGATA        CGATC        CGATG        CGATT        CGATN        CGANA        CGANC 
0.0000000000 0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0001563454 0.0001563454 
       CGANG        CGANT        CGANN        CGCAA        CGCAC        CGCAG        CGCAT        CGCAN        CGCCA        CGCCC        CGCCG        CGCCT        CGCCN 
0.0001563454 0.0001563454 0.0001563454 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCGA        CGCGC        CGCGG        CGCGT        CGCGN        CGCTA        CGCTC        CGCTG        CGCTT        CGCTN        CGCNA        CGCNC        CGCNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCNT        CGCNN        CGGAA        CGGAC        CGGAG        CGGAT        CGGAN        CGGCA        CGGCC        CGGCG        CGGCT        CGGCN        CGGGA 
0.0000000000 0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 0.0010441708 0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 
       CGGGC        CGGGG        CGGGT        CGGGN        CGGTA        CGGTC        CGGTG        CGGTT        CGGTN        CGGNA        CGGNC        CGGNG        CGGNT 
0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0048718941 0.0000000000 0.0005768219 0.0024290355 
       CGGNN        CGTAA        CGTAC        CGTAG        CGTAT        CGTAN        CGTCA        CGTCC        CGTCG        CGTCT        CGTCN        CGTGA        CGTGC 
0.0019694379 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0016522853 0.0016522853 
       CGTGG        CGTGT        CGTGN        CGTTA        CGTTC        CGTTG        CGTTT        CGTTN        CGTNA        CGTNC        CGTNG        CGTNT        CGTNN 
0.0016522853 0.0016522853 0.0016522853 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0004130713 0.0004130713 0.0004130713 0.0004130713 0.0004130713 
       CGNAA        CGNAC        CGNAG        CGNAT        CGNAN        CGNCA        CGNCC        CGNCG        CGNCT        CGNCN        CGNGA        CGNGC        CGNGG 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       CGNGT        CGNGN        CGNTA        CGNTC        CGNTG        CGNTT        CGNTN        CGNNA        CGNNC        CGNNG        CGNNT        CGNNN 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
 [ reached getOption("max.print") -- omitted 2125 entries ]

```


```R

# View the number of S mutations used for estimating mutabilities
model@mutability@numMutS
```


```
[1] 723

```



See also
-------------------

See [TargetingModel](TargetingModel-class.md) for the return object. 
See [plotMutability](plotMutability.md) plotting a mutability model.
See [createSubstitutionMatrix](createSubstitutionMatrix.md), [extendSubstitutionMatrix](extendSubstitutionMatrix.md), 
[createMutabilityMatrix](createMutabilityMatrix.md), [extendMutabilityMatrix](extendMutabilityMatrix.md) and 
[createTargetingMatrix](createTargetingMatrix.md) for component steps in building a model.






