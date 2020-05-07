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
       AAAAA        AAAAC        AAAAG        AAAAT        AAAAN        AAACA        AAACC        AAACG        AAACT        AAACN        AAAGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAAGC        AAAGG        AAAGT        AAAGN        AAATA        AAATC        AAATG        AAATT        AAATN        AAANA        AAANC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAANG        AAANT        AAANN        AACAA        AACAC        AACAG        AACAT        AACAN        AACCA        AACCC        AACCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACCT        AACCN        AACGA        AACGC        AACGG        AACGT        AACGN        AACTA        AACTC        AACTG        AACTT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACTN        AACNA        AACNC        AACNG        AACNT        AACNN        AAGAA        AAGAC        AAGAG        AAGAT        AAGAN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0027781351 0.0000000000 0.0023072876 0.0000000000 0.0012713557 
       AAGCA        AAGCC        AAGCG        AAGCT        AAGCN        AAGGA        AAGGC        AAGGG        AAGGT        AAGGN        AAGTA 
0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0002851050 0.0002066547 0.0148040595 
       AAGTC        AAGTG        AAGTT        AAGTN        AAGNA        AAGNC        AAGNG        AAGNT        AAGNN        AATAA        AATAC 
0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0051708243 0.0000000000 0.0005768219 0.0021539494 0.0019753989 0.0000000000 0.0000000000 
       AATAG        AATAT        AATAN        AATCA        AATCC        AATCG        AATCT        AATCN        AATGA        AATGC        AATGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AATGT        AATGN        AATTA        AATTC        AATTG        AATTT        AATTN        AATNA        AATNC        AATNG        AATNT 
0.0000000000 0.0000000000 0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0001450970 0.0001450970 0.0001450970 0.0001450970 
       AATNN        AANAA        AANAC        AANAG        AANAT        AANAN        AANCA        AANCC        AANCG        AANCT        AANCN 
0.0001450970           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       AANGA        AANGC        AANGG        AANGT        AANGN        AANTA        AANTC        AANTG        AANTT        AANTN        AANNA 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       AANNC        AANNG        AANNT        AANNN        ACAAA        ACAAC        ACAAG        ACAAT        ACAAN        ACACA        ACACC 
          NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACACG        ACACT        ACACN        ACAGA        ACAGC        ACAGG        ACAGT        ACAGN        ACATA        ACATC        ACATG 
0.0000000000 0.0000000000 0.0000000000 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0000000000 0.0000000000 0.0000000000 
       ACATT        ACATN        ACANA        ACANC        ACANG        ACANT        ACANN        ACCAA        ACCAC        ACCAG        ACCAT 
0.0000000000 0.0000000000 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0043115900 0.0043115900 0.0043115900 0.0012186597 
       ACCAN        ACCCA        ACCCC        ACCCG        ACCCT        ACCCN        ACCGA        ACCGC        ACCGG        ACCGT        ACCGN 
0.0035383574 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 
       ACCTA        ACCTC        ACCTG        ACCTT        ACCTN        ACCNA        ACCNC        ACCNG        ACCNT        ACCNN        ACGAA 
0.0012186597 0.0012186597 0.0012186597 0.0043115900 0.0019918923 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0009606561 
       ACGAC        ACGAG        ACGAT        ACGAN        ACGCA        ACGCC        ACGCG        ACGCT        ACGCN        ACGGA        ACGGC 
0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 
       ACGGG        ACGGT        ACGGN        ACGTA        ACGTC        ACGTG        ACGTT        ACGTN        ACGNA        ACGNC        ACGNG 
0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0047164545 0.0000000000 0.0005768219 
       ACGNT        ACGNN        ACTAA        ACTAC        ACTAG        ACTAT        ACTAN        ACTCA        ACTCC        ACTCG        ACTCT 
0.0024290355 0.0019305780 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACTCN        ACTGA        ACTGC        ACTGG        ACTGT        ACTGN        ACTTA        ACTTC        ACTTG        ACTTT        ACTTN 
0.0000000000 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACTNA        ACTNC        ACTNG        ACTNT        ACTNN        ACNAA        ACNAC        ACNAG        ACNAT        ACNAN        ACNCA 
0.0005389060 0.0005389060 0.0005389060 0.0005389060 0.0005389060           NA           NA           NA           NA           NA           NA 
       ACNCC        ACNCG        ACNCT        ACNCN        ACNGA        ACNGC        ACNGG        ACNGT        ACNGN        ACNTA        ACNTC 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ACNTG        ACNTT        ACNTN        ACNNA        ACNNC        ACNNG        ACNNT        ACNNN        AGAAA        AGAAC        AGAAG 
          NA           NA           NA           NA           NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 
       AGAAT        AGAAN        AGACA        AGACC        AGACG        AGACT        AGACN        AGAGA        AGAGC        AGAGG        AGAGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0006253815 0.0006253815 0.0006253815 0.0006253815 
       AGAGN        AGATA        AGATC        AGATG        AGATT        AGATN        AGANA        AGANC        AGANG        AGANT        AGANN 
0.0006253815 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0001563454 0.0001563454 0.0001563454 0.0001563454 0.0001563454 
       AGCAA        AGCAC        AGCAG        AGCAT        AGCAN        AGCCA        AGCCC        AGCCG        AGCCT        AGCCN        AGCGA 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCGC        AGCGG        AGCGT        AGCGN        AGCTA        AGCTC        AGCTG        AGCTT        AGCTN        AGCNA        AGCNC 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCNG        AGCNT        AGCNN        AGGAA        AGGAC        AGGAG        AGGAT        AGGAN        AGGCA        AGGCC        AGGCG 
0.0073813816 0.0073813816 0.0073813816 0.0009606561 0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0025595886 0.0000000000 0.0000000000 
       AGGCT        AGGCN        AGGGA        AGGGC        AGGGG        AGGGT        AGGGN        AGGTA        AGGTC        AGGTG        AGGTT 
0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 
       AGGTN        AGGNA        AGGNC        AGGNG        AGGNT        AGGNN        AGTAA        AGTAC        AGTAG        AGTAT        AGTAN 
0.0037010149 0.0047164545 0.0000000000 0.0005768219 0.0024290355 0.0019305780 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AGTCA        AGTCC        AGTCG        AGTCT        AGTCN        AGTGA        AGTGC        AGTGG        AGTGT        AGTGN        AGTTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0000000000 
       AGTTC        AGTTG        AGTTT        AGTTN        AGTNA        AGTNC        AGTNG        AGTNT        AGTNN        AGNAA        AGNAC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0004130713 0.0004130713 0.0004130713 0.0004130713 0.0004130713           NA           NA 
       AGNAG        AGNAT        AGNAN        AGNCA        AGNCC        AGNCG        AGNCT        AGNCN        AGNGA        AGNGC        AGNGG 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       AGNGT        AGNGN        AGNTA        AGNTC        AGNTG        AGNTT        AGNTN        AGNNA        AGNNC        AGNNG        AGNNT 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       AGNNN        ATAAA        ATAAC        ATAAG        ATAAT        ATAAN        ATACA        ATACC        ATACG        ATACT        ATACN 
          NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATAGA        ATAGC        ATAGG        ATAGT        ATAGN        ATATA        ATATC        ATATG        ATATT        ATATN        ATANA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATANC        ATANG        ATANT        ATANN        ATCAA        ATCAC        ATCAG        ATCAT        ATCAN        ATCCA        ATCCC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCCG        ATCCT        ATCCN        ATCGA        ATCGC        ATCGG        ATCGT        ATCGN        ATCTA        ATCTC        ATCTG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCTT        ATCTN        ATCNA        ATCNC        ATCNG        ATCNT        ATCNN        ATGAA        ATGAC        ATGAG        ATGAT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 
       ATGAN        ATGCA        ATGCC        ATGCG        ATGCT        ATGCN        ATGGA        ATGGC        ATGGG        ATGGT        ATGGN 
0.0010441708 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0024857935 0.0007568269 
       ATGTA        ATGTC        ATGTG        ATGTT        ATGTN        ATGNA        ATGNC        ATGNG        ATGNT        ATGNN        ATTAA 
0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0049436394 0.0000000000 0.0005768219 0.0027041215 0.0020561457 0.0000000000 
       ATTAC        ATTAG        ATTAT        ATTAN        ATTCA        ATTCC        ATTCG        ATTCT        ATTCN        ATTGA        ATTGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTGG        ATTGT        ATTGN        ATTTA        ATTTC        ATTTG        ATTTT        ATTTN        ATTNA        ATTNC        ATTNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTNT        ATTNN        ATNAA        ATNAC        ATNAG        ATNAT        ATNAN        ATNCA        ATNCC        ATNCG        ATNCT 
0.0000000000 0.0000000000           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ATNCN        ATNGA        ATNGC        ATNGG        ATNGT        ATNGN        ATNTA        ATNTC        ATNTG        ATNTT        ATNTN 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ATNNA        ATNNC        ATNNG        ATNNT        ATNNN        ANAAA        ANAAC        ANAAG        ANAAT        ANAAN        ANACA 
          NA           NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANACC        ANACG        ANACT        ANACN        ANAGA        ANAGC        ANAGG        ANAGT        ANAGN        ANATA        ANATC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0010554064 0.0010554064 0.0010554064 0.0010554064 0.0010554064 0.0000000000 0.0000000000 
       ANATG        ANATT        ANATN        ANANA        ANANC        ANANG        ANANT        ANANN        ANCAA        ANCAC        ANCAG 
0.0000000000 0.0000000000 0.0000000000 0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0029232429 0.0029232429 0.0029232429 
       ANCAT        ANCAN        ANCCA        ANCCC        ANCCG        ANCCT        ANCCN        ANCGA        ANCGC        ANCGG        ANCGT 
0.0021500103 0.0027299347 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 
       ANCGN        ANCTA        ANCTC        ANCTG        ANCTT        ANCTN        ANCNA        ANCNC        ANCNG        ANCNT        ANCNN 
0.0025366266 0.0021500103 0.0021500103 0.0021500103 0.0029232429 0.0023433185 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 
       ANGAA        ANGAC        ANGAG        ANGAT        ANGAN        ANGCA        ANGCC        ANGCG        ANGCT        ANGCN        ANGGA 
0.0016422107 0.0000000000 0.0023072876 0.0000000000 0.0009873746 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 
       ANGGC        ANGGG        ANGGT        ANGGN        ANGTA        ANGTC        ANGTG        ANGTT        ANGTN        ANGNA        ANGNC 
0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0048868432 0.0000000000 
       ANGNG        ANGNT        ANGNN        ANTAA        ANTAC        ANTAG        ANTAT        ANTAN        ANTCA        ANTCC        ANTCG 
0.0005768219 0.0024290355 0.0019731751 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANTCT        ANTCN        ANTGA        ANTGC        ANTGG        ANTGT        ANTGN        ANTTA        ANTTC        ANTTG        ANTTT 
0.0000000000 0.0000000000 0.0009519773 0.0009519773 0.0009519773 0.0009519773 0.0009519773 0.0001450970 0.0001450970 0.0001450970 0.0001450970 
       ANTTN        ANTNA        ANTNC        ANTNG        ANTNT        ANTNN        ANNAA        ANNAC        ANNAG        ANNAT        ANNAN 
0.0001450970 0.0002742686 0.0002742686 0.0002742686 0.0002742686 0.0002742686           NA           NA           NA           NA           NA 
       ANNCA        ANNCC        ANNCG        ANNCT        ANNCN        ANNGA        ANNGC        ANNGG        ANNGT        ANNGN        ANNTA 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       ANNTC        ANNTG        ANNTT        ANNTN        ANNNA        ANNNC        ANNNG        ANNNT        ANNNN        CAAAA        CAAAC 
          NA           NA           NA           NA           NA           NA           NA           NA           NA 0.0000000000 0.0000000000 
       CAAAG        CAAAT        CAAAN        CAACA        CAACC        CAACG        CAACT        CAACN        CAAGA        CAAGC        CAAGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAAGT        CAAGN        CAATA        CAATC        CAATG        CAATT        CAATN        CAANA        CAANC        CAANG        CAANT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAANN        CACAA        CACAC        CACAG        CACAT        CACAN        CACCA        CACCC        CACCG        CACCT        CACCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACGA        CACGC        CACGG        CACGT        CACGN        CACTA        CACTC        CACTG        CACTT        CACTN        CACNA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACNC        CACNG        CACNT        CACNN        CAGAA        CAGAC        CAGAG        CAGAT        CAGAN        CAGCA        CAGCC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0009606561 0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0022726072 0.0000000000 
       CAGCG        CAGCT        CAGCN        CAGGA        CAGGC        CAGGG        CAGGT        CAGGN        CAGTA        CAGTC        CAGTG 
0.0000000000 0.0083306927 0.0026508250 0.0005415140 0.0000000000 0.0000000000 0.0024857935 0.0007568269 0.0148040595 0.0000000000 0.0000000000 
       CAGTT        CAGTN        CAGNA        CAGNC        CAGNG        CAGNT        CAGNN        CATAA        CATAC        CATAG        CATAT 
0.0000000000 0.0037010149 0.0046447092 0.0000000000 0.0005768219 0.0027041215 0.0019814132 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CATAN        CATCA        CATCC        CATCG        CATCT        CATCN        CATGA        CATGC        CATGG        CATGT        CATGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CATTA        CATTC        CATTG        CATTT        CATTN        CATNA        CATNC        CATNG        CATNT        CATNN        CANAA 
0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0001450970           NA 
       CANAC        CANAG        CANAT        CANAN        CANCA        CANCC        CANCG        CANCT        CANCN        CANGA        CANGC 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       CANGG        CANGT        CANGN        CANTA        CANTC        CANTG        CANTT        CANTN        CANNA        CANNC        CANNG 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       CANNT        CANNN        CCAAA        CCAAC        CCAAG        CCAAT        CCAAN        CCACA        CCACC        CCACG        CCACT 
          NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCACN        CCAGA        CCAGC        CCAGG        CCAGT        CCAGN        CCATA        CCATC        CCATG        CCATT        CCATN 
0.0000000000 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCANA        CCANC        CCANG        CCANT        CCANN        CCCAA        CCCAC        CCCAG        CCCAT        CCCAN        CCCCA 
0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCCC        CCCCG        CCCCT        CCCCN        CCCGA        CCCGC        CCCGG        CCCGT        CCCGN        CCCTA        CCCTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCTG        CCCTT        CCCTN        CCCNA        CCCNC        CCCNG        CCCNT        CCCNN        CCGAA        CCGAC        CCGAG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0018693956 0.0000000000 0.0023072876 
       CCGAT        CCGAN        CCGCA        CCGCC        CCGCG        CCGCT        CCGCN        CCGGA        CCGGC        CCGGG        CCGGT 
0.0000000000 0.0010441708 0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 0.0000000000 0.0000000000 0.0013854492 
       CCGGN        CCGTA        CCGTC        CCGTG        CCGTT        CCGTN        CCGNA        CCGNC        CCGNG        CCGNT        CCGNN 
0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0048718941 0.0000000000 0.0005768219 0.0024290355 0.0019694379 
       CCTAA        CCTAC        CCTAG        CCTAT        CCTAN        CCTCA        CCTCC        CCTCG        CCTCT        CCTCN        CCTGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0021556238 
       CCTGC        CCTGG        CCTGT        CCTGN        CCTTA        CCTTC        CCTTG        CCTTT        CCTTN        CCTNA        CCTNC 
0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005389060 0.0005389060 
       CCTNG        CCTNT        CCTNN        CCNAA        CCNAC        CCNAG        CCNAT        CCNAN        CCNCA        CCNCC        CCNCG 
0.0005389060 0.0005389060 0.0005389060           NA           NA           NA           NA           NA           NA           NA           NA 
       CCNCT        CCNCN        CCNGA        CCNGC        CCNGG        CCNGT        CCNGN        CCNTA        CCNTC        CCNTG        CCNTT 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       CCNTN        CCNNA        CCNNC        CCNNG        CCNNT        CCNNN        CGAAA        CGAAC        CGAAG        CGAAT        CGAAN 
          NA           NA           NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGACA        CGACC        CGACG        CGACT        CGACN        CGAGA        CGAGC        CGAGG        CGAGT        CGAGN        CGATA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0000000000 
       CGATC        CGATG        CGATT        CGATN        CGANA        CGANC        CGANG        CGANT        CGANN        CGCAA        CGCAC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0001563454 0.0001563454 0.0001563454 0.0001563454 0.0001563454 0.0000000000 0.0000000000 
       CGCAG        CGCAT        CGCAN        CGCCA        CGCCC        CGCCG        CGCCT        CGCCN        CGCGA        CGCGC        CGCGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCGT        CGCGN        CGCTA        CGCTC        CGCTG        CGCTT        CGCTN        CGCNA        CGCNC        CGCNG        CGCNT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCNN        CGGAA        CGGAC        CGGAG        CGGAT        CGGAN        CGGCA        CGGCC        CGGCG        CGGCT        CGGCN 
0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 0.0010441708 0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 
       CGGGA        CGGGC        CGGGG        CGGGT        CGGGN        CGGTA        CGGTC        CGGTG        CGGTT        CGGTN        CGGNA 
0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0048718941 
       CGGNC        CGGNG        CGGNT        CGGNN        CGTAA        CGTAC        CGTAG        CGTAT        CGTAN        CGTCA        CGTCC 
0.0000000000 0.0005768219 0.0024290355 0.0019694379 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGTCG        CGTCT        CGTCN        CGTGA        CGTGC        CGTGG        CGTGT        CGTGN        CGTTA        CGTTC        CGTTG 
0.0000000000 0.0000000000 0.0000000000 0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0000000000 0.0000000000 0.0000000000 
       CGTTT        CGTTN        CGTNA        CGTNC        CGTNG        CGTNT        CGTNN        CGNAA        CGNAC        CGNAG        CGNAT 
0.0000000000 0.0000000000 0.0004130713 0.0004130713 0.0004130713 0.0004130713 0.0004130713           NA           NA           NA           NA 
       CGNAN        CGNCA        CGNCC        CGNCG        CGNCT        CGNCN        CGNGA        CGNGC        CGNGG        CGNGT        CGNGN 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
       CGNTA        CGNTC        CGNTG        CGNTT        CGNTN        CGNNA        CGNNC        CGNNG        CGNNT        CGNNN 
          NA           NA           NA           NA           NA           NA           NA           NA           NA           NA 
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






