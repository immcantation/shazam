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
       AAAAA        AAAAC        AAAAG        AAAAT        AAAAN        AAACA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAACC        AAACG        AAACT        AAACN        AAAGA        AAAGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAAGG        AAAGT        AAAGN        AAATA        AAATC        AAATG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAATT        AAATN        AAANA        AAANC        AAANG        AAANT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAANN        AACAA        AACAC        AACAG        AACAT        AACAN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACCA        AACCC        AACCG        AACCT        AACCN        AACGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACGC        AACGG        AACGT        AACGN        AACTA        AACTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACTG        AACTT        AACTN        AACNA        AACNC        AACNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACNT        AACNN        AAGAA        AAGAC        AAGAG        AAGAT 
0.0000000000 0.0000000000 0.0027781351 0.0000000000 0.0023072876 0.0000000000 
       AAGAN        AAGCA        AAGCC        AAGCG        AAGCT        AAGCN 
0.0012713557 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 
       AAGGA        AAGGC        AAGGG        AAGGT        AAGGN        AAGTA 
0.0005415140 0.0000000000 0.0000000000 0.0002851050 0.0002066547 0.0148040595 
       AAGTC        AAGTG        AAGTT        AAGTN        AAGNA        AAGNC 
0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0051708243 0.0000000000 
       AAGNG        AAGNT        AAGNN        AATAA        AATAC        AATAG 
0.0005768219 0.0021539494 0.0019753989 0.0000000000 0.0000000000 0.0000000000 
       AATAT        AATAN        AATCA        AATCC        AATCG        AATCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AATCN        AATGA        AATGC        AATGG        AATGT        AATGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AATTA        AATTC        AATTG        AATTT        AATTN        AATNA 
0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0001450970 
       AATNC        AATNG        AATNT        AATNN        AANAA        AANAC 
0.0001450970 0.0001450970 0.0001450970 0.0001450970           NA           NA 
       AANAG        AANAT        AANAN        AANCA        AANCC        AANCG 
          NA           NA           NA           NA           NA           NA 
       AANCT        AANCN        AANGA        AANGC        AANGG        AANGT 
          NA           NA           NA           NA           NA           NA 
       AANGN        AANTA        AANTC        AANTG        AANTT        AANTN 
          NA           NA           NA           NA           NA           NA 
       AANNA        AANNC        AANNG        AANNT        AANNN        ACAAA 
          NA           NA           NA           NA           NA 0.0000000000 
       ACAAC        ACAAG        ACAAT        ACAAN        ACACA        ACACC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACACG        ACACT        ACACN        ACAGA        ACAGC        ACAGG 
0.0000000000 0.0000000000 0.0000000000 0.0035962442 0.0035962442 0.0035962442 
       ACAGT        ACAGN        ACATA        ACATC        ACATG        ACATT 
0.0035962442 0.0035962442 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACATN        ACANA        ACANC        ACANG        ACANT        ACANN 
0.0000000000 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0008990611 
       ACCAA        ACCAC        ACCAG        ACCAT        ACCAN        ACCCA 
0.0043115900 0.0043115900 0.0043115900 0.0012186597 0.0035383574 0.0027651248 
       ACCCC        ACCCG        ACCCT        ACCCN        ACCGA        ACCGC 
0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 
       ACCGG        ACCGT        ACCGN        ACCTA        ACCTC        ACCTG 
0.0027651248 0.0027651248 0.0027651248 0.0012186597 0.0012186597 0.0012186597 
       ACCTT        ACCTN        ACCNA        ACCNC        ACCNG        ACCNT 
0.0043115900 0.0019918923 0.0027651248 0.0027651248 0.0027651248 0.0027651248 
       ACCNN        ACGAA        ACGAC        ACGAG        ACGAT        ACGAN 
0.0027651248 0.0009606561 0.0000000000 0.0023072876 0.0000000000 0.0008169859 
       ACGCA        ACGCC        ACGCG        ACGCT        ACGCN        ACGGA 
0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 
       ACGGC        ACGGG        ACGGT        ACGGN        ACGTA        ACGTC 
0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 
       ACGTG        ACGTT        ACGTN        ACGNA        ACGNC        ACGNG 
0.0000000000 0.0000000000 0.0037010149 0.0047164545 0.0000000000 0.0005768219 
       ACGNT        ACGNN        ACTAA        ACTAC        ACTAG        ACTAT 
0.0024290355 0.0019305780 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACTAN        ACTCA        ACTCC        ACTCG        ACTCT        ACTCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACTGA        ACTGC        ACTGG        ACTGT        ACTGN        ACTTA 
0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0000000000 
       ACTTC        ACTTG        ACTTT        ACTTN        ACTNA        ACTNC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005389060 0.0005389060 
       ACTNG        ACTNT        ACTNN        ACNAA        ACNAC        ACNAG 
0.0005389060 0.0005389060 0.0005389060           NA           NA           NA 
       ACNAT        ACNAN        ACNCA        ACNCC        ACNCG        ACNCT 
          NA           NA           NA           NA           NA           NA 
       ACNCN        ACNGA        ACNGC        ACNGG        ACNGT        ACNGN 
          NA           NA           NA           NA           NA           NA 
       ACNTA        ACNTC        ACNTG        ACNTT        ACNTN        ACNNA 
          NA           NA           NA           NA           NA           NA 
       ACNNC        ACNNG        ACNNT        ACNNN        AGAAA        AGAAC 
          NA           NA           NA           NA 0.0000000000 0.0000000000 
       AGAAG        AGAAT        AGAAN        AGACA        AGACC        AGACG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AGACT        AGACN        AGAGA        AGAGC        AGAGG        AGAGT 
0.0000000000 0.0000000000 0.0006253815 0.0006253815 0.0006253815 0.0006253815 
       AGAGN        AGATA        AGATC        AGATG        AGATT        AGATN 
0.0006253815 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AGANA        AGANC        AGANG        AGANT        AGANN        AGCAA 
0.0001563454 0.0001563454 0.0001563454 0.0001563454 0.0001563454 0.0073813816 
       AGCAC        AGCAG        AGCAT        AGCAN        AGCCA        AGCCC 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCCG        AGCCT        AGCCN        AGCGA        AGCGC        AGCGG 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCGT        AGCGN        AGCTA        AGCTC        AGCTG        AGCTT 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCTN        AGCNA        AGCNC        AGCNG        AGCNT        AGCNN 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGGAA        AGGAC        AGGAG        AGGAT        AGGAN        AGGCA 
0.0009606561 0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0025595886 
       AGGCC        AGGCG        AGGCT        AGGCN        AGGGA        AGGGC 
0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 
       AGGGG        AGGGT        AGGGN        AGGTA        AGGTC        AGGTG 
0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 
       AGGTT        AGGTN        AGGNA        AGGNC        AGGNG        AGGNT 
0.0000000000 0.0037010149 0.0047164545 0.0000000000 0.0005768219 0.0024290355 
       AGGNN        AGTAA        AGTAC        AGTAG        AGTAT        AGTAN 
0.0019305780 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AGTCA        AGTCC        AGTCG        AGTCT        AGTCN        AGTGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0016522853 
       AGTGC        AGTGG        AGTGT        AGTGN        AGTTA        AGTTC 
0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0000000000 0.0000000000 
       AGTTG        AGTTT        AGTTN        AGTNA        AGTNC        AGTNG 
0.0000000000 0.0000000000 0.0000000000 0.0004130713 0.0004130713 0.0004130713 
       AGTNT        AGTNN        AGNAA        AGNAC        AGNAG        AGNAT 
0.0004130713 0.0004130713           NA           NA           NA           NA 
       AGNAN        AGNCA        AGNCC        AGNCG        AGNCT        AGNCN 
          NA           NA           NA           NA           NA           NA 
       AGNGA        AGNGC        AGNGG        AGNGT        AGNGN        AGNTA 
          NA           NA           NA           NA           NA           NA 
       AGNTC        AGNTG        AGNTT        AGNTN        AGNNA        AGNNC 
          NA           NA           NA           NA           NA           NA 
       AGNNG        AGNNT        AGNNN        ATAAA        ATAAC        ATAAG 
          NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 
       ATAAT        ATAAN        ATACA        ATACC        ATACG        ATACT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATACN        ATAGA        ATAGC        ATAGG        ATAGT        ATAGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATATA        ATATC        ATATG        ATATT        ATATN        ATANA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATANC        ATANG        ATANT        ATANN        ATCAA        ATCAC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCAG        ATCAT        ATCAN        ATCCA        ATCCC        ATCCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCCT        ATCCN        ATCGA        ATCGC        ATCGG        ATCGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCGN        ATCTA        ATCTC        ATCTG        ATCTT        ATCTN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCNA        ATCNC        ATCNG        ATCNT        ATCNN        ATGAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0018693956 
       ATGAC        ATGAG        ATGAT        ATGAN        ATGCA        ATGCC 
0.0000000000 0.0023072876 0.0000000000 0.0010441708 0.0025595886 0.0000000000 
       ATGCG        ATGCT        ATGCN        ATGGA        ATGGC        ATGGG 
0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 
       ATGGT        ATGGN        ATGTA        ATGTC        ATGTG        ATGTT 
0.0024857935 0.0007568269 0.0148040595 0.0000000000 0.0000000000 0.0000000000 
       ATGTN        ATGNA        ATGNC        ATGNG        ATGNT        ATGNN 
0.0037010149 0.0049436394 0.0000000000 0.0005768219 0.0027041215 0.0020561457 
       ATTAA        ATTAC        ATTAG        ATTAT        ATTAN        ATTCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTCC        ATTCG        ATTCT        ATTCN        ATTGA        ATTGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTGG        ATTGT        ATTGN        ATTTA        ATTTC        ATTTG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTTT        ATTTN        ATTNA        ATTNC        ATTNG        ATTNT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTNN        ATNAA        ATNAC        ATNAG        ATNAT        ATNAN 
0.0000000000           NA           NA           NA           NA           NA 
       ATNCA        ATNCC        ATNCG        ATNCT        ATNCN        ATNGA 
          NA           NA           NA           NA           NA           NA 
       ATNGC        ATNGG        ATNGT        ATNGN        ATNTA        ATNTC 
          NA           NA           NA           NA           NA           NA 
       ATNTG        ATNTT        ATNTN        ATNNA        ATNNC        ATNNG 
          NA           NA           NA           NA           NA           NA 
       ATNNT        ATNNN        ANAAA        ANAAC        ANAAG        ANAAT 
          NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANAAN        ANACA        ANACC        ANACG        ANACT        ANACN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANAGA        ANAGC        ANAGG        ANAGT        ANAGN        ANATA 
0.0010554064 0.0010554064 0.0010554064 0.0010554064 0.0010554064 0.0000000000 
       ANATC        ANATG        ANATT        ANATN        ANANA        ANANC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0002638516 0.0002638516 
       ANANG        ANANT        ANANN        ANCAA        ANCAC        ANCAG 
0.0002638516 0.0002638516 0.0002638516 0.0029232429 0.0029232429 0.0029232429 
       ANCAT        ANCAN        ANCCA        ANCCC        ANCCG        ANCCT 
0.0021500103 0.0027299347 0.0025366266 0.0025366266 0.0025366266 0.0025366266 
       ANCCN        ANCGA        ANCGC        ANCGG        ANCGT        ANCGN 
0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 
       ANCTA        ANCTC        ANCTG        ANCTT        ANCTN        ANCNA 
0.0021500103 0.0021500103 0.0021500103 0.0029232429 0.0023433185 0.0025366266 
       ANCNC        ANCNG        ANCNT        ANCNN        ANGAA        ANGAC 
0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0016422107 0.0000000000 
       ANGAG        ANGAT        ANGAN        ANGCA        ANGCC        ANGCG 
0.0023072876 0.0000000000 0.0009873746 0.0025595886 0.0000000000 0.0000000000 
       ANGCT        ANGCN        ANGGA        ANGGC        ANGGG        ANGGT 
0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 
       ANGGN        ANGTA        ANGTC        ANGTG        ANGTT        ANGTN 
0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 
       ANGNA        ANGNC        ANGNG        ANGNT        ANGNN        ANTAA 
0.0048868432 0.0000000000 0.0005768219 0.0024290355 0.0019731751 0.0000000000 
       ANTAC        ANTAG        ANTAT        ANTAN        ANTCA        ANTCC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANTCG        ANTCT        ANTCN        ANTGA        ANTGC        ANTGG 
0.0000000000 0.0000000000 0.0000000000 0.0009519773 0.0009519773 0.0009519773 
       ANTGT        ANTGN        ANTTA        ANTTC        ANTTG        ANTTT 
0.0009519773 0.0009519773 0.0001450970 0.0001450970 0.0001450970 0.0001450970 
       ANTTN        ANTNA        ANTNC        ANTNG        ANTNT        ANTNN 
0.0001450970 0.0002742686 0.0002742686 0.0002742686 0.0002742686 0.0002742686 
       ANNAA        ANNAC        ANNAG        ANNAT        ANNAN        ANNCA 
          NA           NA           NA           NA           NA           NA 
       ANNCC        ANNCG        ANNCT        ANNCN        ANNGA        ANNGC 
          NA           NA           NA           NA           NA           NA 
       ANNGG        ANNGT        ANNGN        ANNTA        ANNTC        ANNTG 
          NA           NA           NA           NA           NA           NA 
       ANNTT        ANNTN        ANNNA        ANNNC        ANNNG        ANNNT 
          NA           NA           NA           NA           NA           NA 
       ANNNN        CAAAA        CAAAC        CAAAG        CAAAT        CAAAN 
          NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAACA        CAACC        CAACG        CAACT        CAACN        CAAGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAAGC        CAAGG        CAAGT        CAAGN        CAATA        CAATC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAATG        CAATT        CAATN        CAANA        CAANC        CAANG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAANT        CAANN        CACAA        CACAC        CACAG        CACAT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACAN        CACCA        CACCC        CACCG        CACCT        CACCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACGA        CACGC        CACGG        CACGT        CACGN        CACTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACTC        CACTG        CACTT        CACTN        CACNA        CACNC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACNG        CACNT        CACNN        CAGAA        CAGAC        CAGAG 
0.0000000000 0.0000000000 0.0000000000 0.0009606561 0.0000000000 0.0023072876 
       CAGAT        CAGAN        CAGCA        CAGCC        CAGCG        CAGCT 
0.0000000000 0.0008169859 0.0022726072 0.0000000000 0.0000000000 0.0083306927 
       CAGCN        CAGGA        CAGGC        CAGGG        CAGGT        CAGGN 
0.0026508250 0.0005415140 0.0000000000 0.0000000000 0.0024857935 0.0007568269 
       CAGTA        CAGTC        CAGTG        CAGTT        CAGTN        CAGNA 
0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0046447092 
       CAGNC        CAGNG        CAGNT        CAGNN        CATAA        CATAC 
0.0000000000 0.0005768219 0.0027041215 0.0019814132 0.0000000000 0.0000000000 
       CATAG        CATAT        CATAN        CATCA        CATCC        CATCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CATCT        CATCN        CATGA        CATGC        CATGG        CATGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CATGN        CATTA        CATTC        CATTG        CATTT        CATTN 
0.0000000000 0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0005803880 
       CATNA        CATNC        CATNG        CATNT        CATNN        CANAA 
0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0001450970           NA 
       CANAC        CANAG        CANAT        CANAN        CANCA        CANCC 
          NA           NA           NA           NA           NA           NA 
       CANCG        CANCT        CANCN        CANGA        CANGC        CANGG 
          NA           NA           NA           NA           NA           NA 
       CANGT        CANGN        CANTA        CANTC        CANTG        CANTT 
          NA           NA           NA           NA           NA           NA 
       CANTN        CANNA        CANNC        CANNG        CANNT        CANNN 
          NA           NA           NA           NA           NA           NA 
       CCAAA        CCAAC        CCAAG        CCAAT        CCAAN        CCACA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCACC        CCACG        CCACT        CCACN        CCAGA        CCAGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0035962442 0.0035962442 
       CCAGG        CCAGT        CCAGN        CCATA        CCATC        CCATG 
0.0035962442 0.0035962442 0.0035962442 0.0000000000 0.0000000000 0.0000000000 
       CCATT        CCATN        CCANA        CCANC        CCANG        CCANT 
0.0000000000 0.0000000000 0.0008990611 0.0008990611 0.0008990611 0.0008990611 
       CCANN        CCCAA        CCCAC        CCCAG        CCCAT        CCCAN 
0.0008990611 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCCA        CCCCC        CCCCG        CCCCT        CCCCN        CCCGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCGC        CCCGG        CCCGT        CCCGN        CCCTA        CCCTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCTG        CCCTT        CCCTN        CCCNA        CCCNC        CCCNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCNT        CCCNN        CCGAA        CCGAC        CCGAG        CCGAT 
0.0000000000 0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 
       CCGAN        CCGCA        CCGCC        CCGCG        CCGCT        CCGCN 
0.0010441708 0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 
       CCGGA        CCGGC        CCGGG        CCGGT        CCGGN        CCGTA 
0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 
       CCGTC        CCGTG        CCGTT        CCGTN        CCGNA        CCGNC 
0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0048718941 0.0000000000 
       CCGNG        CCGNT        CCGNN        CCTAA        CCTAC        CCTAG 
0.0005768219 0.0024290355 0.0019694379 0.0000000000 0.0000000000 0.0000000000 
       CCTAT        CCTAN        CCTCA        CCTCC        CCTCG        CCTCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCTCN        CCTGA        CCTGC        CCTGG        CCTGT        CCTGN 
0.0000000000 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0021556238 
       CCTTA        CCTTC        CCTTG        CCTTT        CCTTN        CCTNA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005389060 
       CCTNC        CCTNG        CCTNT        CCTNN        CCNAA        CCNAC 
0.0005389060 0.0005389060 0.0005389060 0.0005389060           NA           NA 
       CCNAG        CCNAT        CCNAN        CCNCA        CCNCC        CCNCG 
          NA           NA           NA           NA           NA           NA 
       CCNCT        CCNCN        CCNGA        CCNGC        CCNGG        CCNGT 
          NA           NA           NA           NA           NA           NA 
       CCNGN        CCNTA        CCNTC        CCNTG        CCNTT        CCNTN 
          NA           NA           NA           NA           NA           NA 
       CCNNA        CCNNC        CCNNG        CCNNT        CCNNN        CGAAA 
          NA           NA           NA           NA           NA 0.0000000000 
       CGAAC        CGAAG        CGAAT        CGAAN        CGACA        CGACC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGACG        CGACT        CGACN        CGAGA        CGAGC        CGAGG 
0.0000000000 0.0000000000 0.0000000000 0.0006253815 0.0006253815 0.0006253815 
       CGAGT        CGAGN        CGATA        CGATC        CGATG        CGATT 
0.0006253815 0.0006253815 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGATN        CGANA        CGANC        CGANG        CGANT        CGANN 
0.0000000000 0.0001563454 0.0001563454 0.0001563454 0.0001563454 0.0001563454 
       CGCAA        CGCAC        CGCAG        CGCAT        CGCAN        CGCCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCCC        CGCCG        CGCCT        CGCCN        CGCGA        CGCGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCGG        CGCGT        CGCGN        CGCTA        CGCTC        CGCTG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCTT        CGCTN        CGCNA        CGCNC        CGCNG        CGCNT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCNN        CGGAA        CGGAC        CGGAG        CGGAT        CGGAN 
0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 0.0010441708 
       CGGCA        CGGCC        CGGCG        CGGCT        CGGCN        CGGGA 
0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 
       CGGGC        CGGGG        CGGGT        CGGGN        CGGTA        CGGTC 
0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 
       CGGTG        CGGTT        CGGTN        CGGNA        CGGNC        CGGNG 
0.0000000000 0.0000000000 0.0037010149 0.0048718941 0.0000000000 0.0005768219 
       CGGNT        CGGNN        CGTAA        CGTAC        CGTAG        CGTAT 
0.0024290355 0.0019694379 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGTAN        CGTCA        CGTCC        CGTCG        CGTCT        CGTCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGTGA        CGTGC        CGTGG        CGTGT        CGTGN        CGTTA 
0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0000000000 
       CGTTC        CGTTG        CGTTT        CGTTN        CGTNA        CGTNC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0004130713 0.0004130713 
       CGTNG        CGTNT        CGTNN        CGNAA        CGNAC        CGNAG 
0.0004130713 0.0004130713 0.0004130713           NA           NA           NA 
       CGNAT        CGNAN        CGNCA        CGNCC        CGNCG        CGNCT 
          NA           NA           NA           NA           NA           NA 
       CGNCN        CGNGA        CGNGC        CGNGG        CGNGT        CGNGN 
          NA           NA           NA           NA           NA           NA 
       CGNTA        CGNTC        CGNTG        CGNTT        CGNTN        CGNNA 
          NA           NA           NA           NA           NA           NA 
       CGNNC        CGNNG        CGNNT        CGNNN        CTAAA        CTAAC 
          NA           NA           NA           NA 0.0000000000 0.0000000000 
       CTAAG        CTAAT        CTAAN        CTACA        CTACC        CTACG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTACT        CTACN        CTAGA        CTAGC        CTAGG        CTAGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTAGN        CTATA        CTATC        CTATG        CTATT        CTATN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTANA        CTANC        CTANG        CTANT        CTANN        CTCAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTCAC        CTCAG        CTCAT        CTCAN        CTCCA        CTCCC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTCCG        CTCCT        CTCCN        CTCGA        CTCGC        CTCGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTCGT        CTCGN        CTCTA        CTCTC        CTCTG        CTCTT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTCTN        CTCNA        CTCNC        CTCNG        CTCNT        CTCNN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTGAA        CTGAC        CTGAG        CTGAT        CTGAN        CTGCA 
0.0027781351 0.0000000000 0.0023072876 0.0000000000 0.0012713557 0.0028465699 
       CTGCC        CTGCG        CTGCT        CTGCN        CTGGA        CTGGC 
0.0000000000 0.0000000000 0.0083306927 0.0027943157 0.0005415140 0.0000000000 
       CTGGG        CTGGT        CTGGN        CTGTA        CTGTC        CTGTG 
0.0000000000 0.0002851050 0.0002066547 0.0148040595 0.0000000000 0.0000000000 
       CTGTT        CTGTN        CTGNA        CTGNC        CTGNG        CTGNT 
0.0000000000 0.0037010149 0.0052425696 0.0000000000 0.0005768219 0.0021539494 
       CTGNN        CTTAA        CTTAC        CTTAG        CTTAT        CTTAN 
0.0019933352 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTTCA        CTTCC        CTTCG        CTTCT        CTTCN        CTTGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTTGC        CTTGG        CTTGT        CTTGN        CTTTA        CTTTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTTTG        CTTTT        CTTTN        CTTNA        CTTNC        CTTNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CTTNT        CTTNN        CTNAA        CTNAC        CTNAG        CTNAT 
0.0000000000 0.0000000000           NA           NA           NA           NA 
       CTNAN        CTNCA        CTNCC        CTNCG        CTNCT        CTNCN 
          NA           NA           NA           NA           NA           NA 
       CTNGA        CTNGC        CTNGG        CTNGT        CTNGN        CTNTA 
          NA           NA           NA           NA           NA           NA 
       CTNTC        CTNTG        CTNTT        CTNTN        CTNNA        CTNNC 
          NA           NA           NA           NA           NA           NA 
       CTNNG        CTNNT        CTNNN        CNAAA        CNAAC        CNAAG 
          NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 
       CNAAT        CNAAN        CNACA        CNACC        CNACG        CNACT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CNACN        CNAGA        CNAGC        CNAGG        CNAGT        CNAGN 
0.0000000000 0.0010554064 0.0010554064 0.0010554064 0.0010554064 0.0010554064 
       CNATA        CNATC        CNATG        CNATT        CNATN        CNANA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0002638516 
       CNANC        CNANG        CNANT        CNANN        CNCAA        CNCAC 
0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0000000000 0.0000000000 
       CNCAG        CNCAT        CNCAN        CNCCA        CNCCC        CNCCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CNCCT        CNCCN        CNCGA        CNCGC        CNCGG        CNCGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CNCGN        CNCTA        CNCTC        CNCTG        CNCTT        CNCTN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CNCNA        CNCNC        CNCNG        CNCNT        CNCNN        CNGAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0018693956 
       CNGAC        CNGAG        CNGAT        CNGAN        CNGCA        CNGCC 
0.0000000000 0.0023072876 0.0000000000 0.0010441708 0.0024160979 0.0000000000 
       CNGCG        CNGCT        CNGCN        CNGGA        CNGGC        CNGGG 
0.0000000000 0.0083306927 0.0026866977 0.0005415140 0.0000000000 0.0000000000 
       CNGGT        CNGGN        CNGTA        CNGTC        CNGTG        CNGTT 
0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 
       CNGTN        CNGNA        CNGNC        CNGNG        CNGNT        CNGNN 
0.0037010149 0.0049077668 0.0000000000 0.0005768219 0.0024290355 0.0019784060 
       CNTAA        CNTAC        CNTAG        CNTAT        CNTAN        CNTCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CNTCC        CNTCG        CNTCT        CNTCN        CNTGA        CNTGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0009519773 0.0009519773 
       CNTGG        CNTGT        CNTGN        CNTTA        CNTTC        CNTTG 
0.0009519773 0.0009519773 0.0009519773 0.0001450970 0.0001450970 0.0001450970 
       CNTTT        CNTTN        CNTNA        CNTNC        CNTNG        CNTNT 
0.0001450970 0.0001450970 0.0002742686 0.0002742686 0.0002742686 0.0002742686 
       CNTNN        CNNAA        CNNAC        CNNAG        CNNAT        CNNAN 
0.0002742686           NA           NA           NA           NA           NA 
       CNNCA        CNNCC        CNNCG        CNNCT        CNNCN        CNNGA 
          NA           NA           NA           NA           NA           NA 
       CNNGC        CNNGG        CNNGT        CNNGN        CNNTA        CNNTC 
          NA           NA           NA           NA           NA           NA 
       CNNTG        CNNTT        CNNTN        CNNNA        CNNNC        CNNNG 
          NA           NA           NA           NA           NA           NA 
       CNNNT        CNNNN        GAAAA        GAAAC        GAAAG        GAAAT 
          NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GAAAN        GAACA        GAACC        GAACG        GAACT        GAACN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GAAGA        GAAGC        GAAGG        GAAGT        GAAGN        GAATA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GAATC        GAATG        GAATT        GAATN        GAANA        GAANC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GAANG        GAANT        GAANN        GACAA        GACAC        GACAG 
0.0000000000 0.0000000000 0.0000000000 0.0005310128 0.0005310128 0.0005310128 
       GACAT        GACAN        GACCA        GACCC        GACCG        GACCT 
0.0005310128 0.0005310128 0.0005310128 0.0005310128 0.0005310128 0.0005310128 
       GACCN        GACGA        GACGC        GACGG        GACGT        GACGN 
0.0005310128 0.0005310128 0.0005310128 0.0005310128 0.0005310128 0.0005310128 
       GACTA        GACTC        GACTG        GACTT        GACTN        GACNA 
0.0005310128 0.0005310128 0.0005310128 0.0005310128 0.0005310128 0.0005310128 
       GACNC        GACNG        GACNT        GACNN        GAGAA        GAGAC 
0.0005310128 0.0005310128 0.0005310128 0.0005310128 0.0018693956 0.0000000000 
       GAGAG        GAGAT        GAGAN        GAGCA        GAGCC        GAGCG 
0.0023072876 0.0000000000 0.0010441708 0.0028465699 0.0000000000 0.0000000000 
       GAGCT        GAGCN        GAGGA        GAGGC        GAGGG        GAGGT 
0.0083306927 0.0027943157 0.0005415140 0.0000000000 0.0000000000 0.0002851050 
       GAGGN        GAGTA        GAGTC        GAGTG        GAGTT        GAGTN 
0.0002066547 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 
       GAGNA        GAGNC        GAGNG        GAGNT        GAGNN        GATAA 
0.0050153848 0.0000000000 0.0005768219 0.0021539494 0.0019365390 0.0000000000 
       GATAC        GATAG        GATAT        GATAN        GATCA        GATCC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GATCG        GATCT        GATCN        GATGA        GATGC        GATGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GATGT        GATGN        GATTA        GATTC        GATTG        GATTT 
0.0000000000 0.0000000000 0.0005803880 0.0005803880 0.0005803880 0.0005803880 
       GATTN        GATNA        GATNC        GATNG        GATNT        GATNN 
0.0005803880 0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0001450970 
       GANAA        GANAC        GANAG        GANAT        GANAN        GANCA 
          NA           NA           NA           NA           NA           NA 
       GANCC        GANCG        GANCT        GANCN        GANGA        GANGC 
          NA           NA           NA           NA           NA           NA 
       GANGG        GANGT        GANGN        GANTA        GANTC        GANTG 
          NA           NA           NA           NA           NA           NA 
       GANTT        GANTN        GANNA        GANNC        GANNG        GANNT 
          NA           NA           NA           NA           NA           NA 
       GANNN        GCAAA        GCAAC        GCAAG        GCAAT        GCAAN 
          NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GCACA        GCACC        GCACG        GCACT        GCACN        GCAGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0035962442 
       GCAGC        GCAGG        GCAGT        GCAGN        GCATA        GCATC 
0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0000000000 0.0000000000 
       GCATG        GCATT        GCATN        GCANA        GCANC        GCANG 
0.0000000000 0.0000000000 0.0000000000 0.0008990611 0.0008990611 0.0008990611 
       GCANT        GCANN        GCCAA        GCCAC        GCCAG        GCCAT 
0.0008990611 0.0008990611 0.0004326306 0.0004326306 0.0004326306 0.0004326306 
       GCCAN        GCCCA        GCCCC        GCCCG        GCCCT        GCCCN 
0.0004326306 0.0004326306 0.0004326306 0.0004326306 0.0004326306 0.0004326306 
       GCCGA        GCCGC        GCCGG        GCCGT        GCCGN        GCCTA 
0.0004326306 0.0004326306 0.0004326306 0.0004326306 0.0004326306 0.0004326306 
       GCCTC        GCCTG        GCCTT        GCCTN        GCCNA        GCCNC 
0.0004326306 0.0004326306 0.0004326306 0.0004326306 0.0004326306 0.0004326306 
       GCCNG        GCCNT        GCCNN        GCGAA        GCGAC        GCGAG 
0.0004326306 0.0004326306 0.0004326306 0.0027781351 0.0000000000 0.0023072876 
       GCGAT        GCGAN        GCGCA        GCGCC        GCGCG        GCGCT 
0.0000000000 0.0012713557 0.0028465699 0.0000000000 0.0000000000 0.0083306927 
       GCGCN        GCGGA        GCGGC        GCGGG        GCGGT        GCGGN 
0.0027943157 0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 
       GCGTA        GCGTC        GCGTG        GCGTT        GCGTN        GCGNA 
0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0052425696 
       GCGNC        GCGNG        GCGNT        GCGNN        GCTAA        GCTAC 
0.0000000000 0.0005768219 0.0024290355 0.0020621068 0.0000000000 0.0000000000 
       GCTAG        GCTAT        GCTAN        GCTCA        GCTCC        GCTCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GCTCT        GCTCN        GCTGA        GCTGC        GCTGG        GCTGT 
0.0000000000 0.0000000000 0.0021556238 0.0021556238 0.0021556238 0.0021556238 
       GCTGN        GCTTA        GCTTC        GCTTG        GCTTT        GCTTN 
0.0021556238 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GCTNA        GCTNC        GCTNG        GCTNT        GCTNN        GCNAA 
0.0005389060 0.0005389060 0.0005389060 0.0005389060 0.0005389060           NA 
       GCNAC        GCNAG        GCNAT        GCNAN        GCNCA        GCNCC 
          NA           NA           NA           NA           NA           NA 
       GCNCG        GCNCT        GCNCN        GCNGA        GCNGC        GCNGG 
          NA           NA           NA           NA           NA           NA 
       GCNGT        GCNGN        GCNTA        GCNTC        GCNTG        GCNTT 
          NA           NA           NA           NA           NA           NA 
       GCNTN        GCNNA        GCNNC        GCNNG        GCNNT        GCNNN 
          NA           NA           NA           NA           NA           NA 
       GGAAA        GGAAC        GGAAG        GGAAT        GGAAN        GGACA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GGACC        GGACG        GGACT        GGACN        GGAGA        GGAGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0006253815 0.0006253815 
       GGAGG        GGAGT        GGAGN        GGATA        GGATC        GGATG 
0.0006253815 0.0006253815 0.0006253815 0.0000000000 0.0000000000 0.0000000000 
       GGATT        GGATN        GGANA        GGANC        GGANG        GGANT 
0.0000000000 0.0000000000 0.0001563454 0.0001563454 0.0001563454 0.0001563454 
       GGANN        GGCAA        GGCAC        GGCAG        GGCAT        GGCAN 
0.0001563454 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GGCCA        GGCCC        GGCCG        GGCCT        GGCCN        GGCGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GGCGC        GGCGG        GGCGT        GGCGN        GGCTA        GGCTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GGCTG        GGCTT        GGCTN        GGCNA        GGCNC        GGCNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GGCNT        GGCNN        GGGAA        GGGAC        GGGAG        GGGAT 
0.0000000000 0.0000000000 0.0027781351 0.0000000000 0.0023072876 0.0000000000 
       GGGAN        GGGCA        GGGCC        GGGCG        GGGCT        GGGCN 
0.0012713557 0.0028465699 0.0000000000 0.0000000000 0.0083306927 0.0027943157 
       GGGGA        GGGGC        GGGGG        GGGGT        GGGGN        GGGTA 
0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 
       GGGTC        GGGTG        GGGTT        GGGTN        GGGNA        GGGNC 
0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0052425696 0.0000000000 
       GGGNG        GGGNT        GGGNN        GGTAA        GGTAC        GGTAG 
0.0005768219 0.0024290355 0.0020621068 0.0000000000 0.0000000000 0.0000000000 
       GGTAT        GGTAN        GGTCA        GGTCC        GGTCG        GGTCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GGTCN        GGTGA        GGTGC        GGTGG        GGTGT        GGTGN 
0.0000000000 0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0016522853 
       GGTTA        GGTTC        GGTTG        GGTTT        GGTTN        GGTNA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0004130713 
       GGTNC        GGTNG        GGTNT        GGTNN        GGNAA        GGNAC 
0.0004130713 0.0004130713 0.0004130713 0.0004130713           NA           NA 
       GGNAG        GGNAT        GGNAN        GGNCA        GGNCC        GGNCG 
          NA           NA           NA           NA           NA           NA 
       GGNCT        GGNCN        GGNGA        GGNGC        GGNGG        GGNGT 
          NA           NA           NA           NA           NA           NA 
       GGNGN        GGNTA        GGNTC        GGNTG        GGNTT        GGNTN 
          NA           NA           NA           NA           NA           NA 
       GGNNA        GGNNC        GGNNG        GGNNT        GGNNN        GTAAA 
          NA           NA           NA           NA           NA 0.0000000000 
       GTAAC        GTAAG        GTAAT        GTAAN        GTACA        GTACC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTACG        GTACT        GTACN        GTAGA        GTAGC        GTAGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTAGT        GTAGN        GTATA        GTATC        GTATG        GTATT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTATN        GTANA        GTANC        GTANG        GTANT        GTANN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTCAA        GTCAC        GTCAG        GTCAT        GTCAN        GTCCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTCCC        GTCCG        GTCCT        GTCCN        GTCGA        GTCGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTCGG        GTCGT        GTCGN        GTCTA        GTCTC        GTCTG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTCTT        GTCTN        GTCNA        GTCNC        GTCNG        GTCNT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTCNN        GTGAA        GTGAC        GTGAG        GTGAT        GTGAN 
0.0000000000 0.0009606561 0.0000000000 0.0023072876 0.0000000000 0.0008169859 
       GTGCA        GTGCC        GTGCG        GTGCT        GTGCN        GTGGA 
0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 
       GTGGC        GTGGG        GTGGT        GTGGN        GTGTA        GTGTC 
0.0000000000 0.0000000000 0.0024857935 0.0007568269 0.0148040595 0.0000000000 
       GTGTG        GTGTT        GTGTN        GTGNA        GTGNC        GTGNG 
0.0000000000 0.0000000000 0.0037010149 0.0046447092 0.0000000000 0.0005768219 
       GTGNT        GTGNN        GTTAA        GTTAC        GTTAG        GTTAT 
0.0027041215 0.0019814132 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTTAN        GTTCA        GTTCC        GTTCG        GTTCT        GTTCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTTGA        GTTGC        GTTGG        GTTGT        GTTGN        GTTTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTTTC        GTTTG        GTTTT        GTTTN        GTTNA        GTTNC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GTTNG        GTTNT        GTTNN        GTNAA        GTNAC        GTNAG 
0.0000000000 0.0000000000 0.0000000000           NA           NA           NA 
       GTNAT        GTNAN        GTNCA        GTNCC        GTNCG        GTNCT 
          NA           NA           NA           NA           NA           NA 
       GTNCN        GTNGA        GTNGC        GTNGG        GTNGT        GTNGN 
          NA           NA           NA           NA           NA           NA 
       GTNTA        GTNTC        GTNTG        GTNTT        GTNTN        GTNNA 
          NA           NA           NA           NA           NA           NA 
       GTNNC        GTNNG        GTNNT        GTNNN        GNAAA        GNAAC 
          NA           NA           NA           NA 0.0000000000 0.0000000000 
       GNAAG        GNAAT        GNAAN        GNACA        GNACC        GNACG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GNACT        GNACN        GNAGA        GNAGC        GNAGG        GNAGT 
0.0000000000 0.0000000000 0.0010554064 0.0010554064 0.0010554064 0.0010554064 
       GNAGN        GNATA        GNATC        GNATG        GNATT        GNATN 
0.0010554064 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GNANA        GNANC        GNANG        GNANT        GNANN        GNCAA 
0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0002409109 
       GNCAC        GNCAG        GNCAT        GNCAN        GNCCA        GNCCC 
0.0002409109 0.0002409109 0.0002409109 0.0002409109 0.0002409109 0.0002409109 
       GNCCG        GNCCT        GNCCN        GNCGA        GNCGC        GNCGG 
0.0002409109 0.0002409109 0.0002409109 0.0002409109 0.0002409109 0.0002409109 
       GNCGT        GNCGN        GNCTA        GNCTC        GNCTG        GNCTT 
0.0002409109 0.0002409109 0.0002409109 0.0002409109 0.0002409109 0.0002409109 
       GNCTN        GNCNA        GNCNC        GNCNG        GNCNT        GNCNN 
0.0002409109 0.0002409109 0.0002409109 0.0002409109 0.0002409109 0.0002409109 
       GNGAA        GNGAC        GNGAG        GNGAT        GNGAN        GNGCA 
0.0020965805 0.0000000000 0.0023072876 0.0000000000 0.0011009670 0.0027030792 
       GNGCC        GNGCG        GNGCT        GNGCN        GNGGA        GNGGC 
0.0000000000 0.0000000000 0.0083306927 0.0027584430 0.0005415140 0.0000000000 
       GNGGG        GNGGT        GNGGN        GNGTA        GNGTC        GNGTG 
0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 
       GNGTT        GNGTN        GNGNA        GNGNC        GNGNG        GNGNT 
0.0000000000 0.0037010149 0.0050363083 0.0000000000 0.0005768219 0.0024290355 
       GNGNN        GNTAA        GNTAC        GNTAG        GNTAT        GNTAN 
0.0020105414 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       GNTCA        GNTCC        GNTCG        GNTCT        GNTCN        GNTGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0009519773 
       GNTGC        GNTGG        GNTGT        GNTGN        GNTTA        GNTTC 
0.0009519773 0.0009519773 0.0009519773 0.0009519773 0.0001450970 0.0001450970 
       GNTTG        GNTTT        GNTTN        GNTNA        GNTNC        GNTNG 
0.0001450970 0.0001450970 0.0001450970 0.0002742686 0.0002742686 0.0002742686 
       GNTNT        GNTNN        GNNAA        GNNAC        GNNAG        GNNAT 
0.0002742686 0.0002742686           NA           NA           NA           NA 
       GNNAN        GNNCA        GNNCC        GNNCG        GNNCT        GNNCN 
          NA           NA           NA           NA           NA           NA 
       GNNGA        GNNGC        GNNGG        GNNGT        GNNGN        GNNTA 
          NA           NA           NA           NA           NA           NA 
       GNNTC        GNNTG        GNNTT        GNNTN        GNNNA        GNNNC 
          NA           NA           NA           NA           NA           NA 
       GNNNG        GNNNT        GNNNN        TAAAA        TAAAC        TAAAG 
          NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 
       TAAAT        TAAAN        TAACA        TAACC        TAACG        TAACT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TAACN        TAAGA        TAAGC        TAAGG        TAAGT        TAAGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TAATA        TAATC        TAATG        TAATT        TAATN        TAANA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TAANC        TAANG        TAANT        TAANN        TACAA        TACAC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0109819402 0.0109819402 
       TACAG        TACAT        TACAN        TACCA        TACCC        TACCG 
0.0109819402 0.0109819402 0.0109819402 0.0109819402 0.0109819402 0.0109819402 
       TACCT        TACCN        TACGA        TACGC        TACGG        TACGT 
0.0109819402 0.0109819402 0.0109819402 0.0109819402 0.0109819402 0.0109819402 
       TACGN        TACTA        TACTC        TACTG        TACTT        TACTN 
0.0109819402 0.0109819402 0.0109819402 0.0109819402 0.0109819402 0.0109819402 
       TACNA        TACNC        TACNG        TACNT        TACNN        TAGAA 
0.0109819402 0.0109819402 0.0109819402 0.0109819402 0.0109819402 0.0009606561 
       TAGAC        TAGAG        TAGAT        TAGAN        TAGCA        TAGCC 
0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0025595886 0.0000000000 
       TAGCG        TAGCT        TAGCN        TAGGA        TAGGC        TAGGG 
0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 
       TAGGT        TAGGN        TAGTA        TAGTC        TAGTG        TAGTT 
0.0002851050 0.0002066547 0.0148040595 0.0000000000 0.0000000000 0.0000000000 
       TAGTN        TAGNA        TAGNC        TAGNG        TAGNT        TAGNN 
0.0037010149 0.0047164545 0.0000000000 0.0005768219 0.0021539494 0.0018618065 
       TATAA        TATAC        TATAG        TATAT        TATAN        TATCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TATCC        TATCG        TATCT        TATCN        TATGA        TATGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TATGG        TATGT        TATGN        TATTA        TATTC        TATTG 
0.0000000000 0.0000000000 0.0000000000 0.0005803880 0.0005803880 0.0005803880 
       TATTT        TATTN        TATNA        TATNC        TATNG        TATNT 
0.0005803880 0.0005803880 0.0001450970 0.0001450970 0.0001450970 0.0001450970 
       TATNN        TANAA        TANAC        TANAG        TANAT        TANAN 
0.0001450970           NA           NA           NA           NA           NA 
       TANCA        TANCC        TANCG        TANCT        TANCN        TANGA 
          NA           NA           NA           NA           NA           NA 
       TANGC        TANGG        TANGT        TANGN        TANTA        TANTC 
          NA           NA           NA           NA           NA           NA 
       TANTG        TANTT        TANTN        TANNA        TANNC        TANNG 
          NA           NA           NA           NA           NA           NA 
       TANNT        TANNN        TCAAA        TCAAC        TCAAG        TCAAT 
          NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TCAAN        TCACA        TCACC        TCACG        TCACT        TCACN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TCAGA        TCAGC        TCAGG        TCAGT        TCAGN        TCATA 
0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0000000000 
       TCATC        TCATG        TCATT        TCATN        TCANA        TCANC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0008990611 0.0008990611 
       TCANG        TCANT        TCANN        TCCAA        TCCAC        TCCAG 
0.0008990611 0.0008990611 0.0008990611 0.0000000000 0.0000000000 0.0000000000 
       TCCAT        TCCAN        TCCCA        TCCCC        TCCCG        TCCCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TCCCN        TCCGA        TCCGC        TCCGG        TCCGT        TCCGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TCCTA        TCCTC        TCCTG        TCCTT        TCCTN        TCCNA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TCCNC        TCCNG        TCCNT        TCCNN        TCGAA        TCGAC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0018693956 0.0000000000 
       TCGAG        TCGAT        TCGAN        TCGCA        TCGCC        TCGCG 
0.0023072876 0.0000000000 0.0010441708 0.0025595886 0.0000000000 0.0000000000 
       TCGCT        TCGCN        TCGGA        TCGGC        TCGGG        TCGGT 
0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 
       TCGGN        TCGTA        TCGTC        TCGTG        TCGTT        TCGTN 
0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 
       TCGNA        TCGNC        TCGNG        TCGNT        TCGNN        TCTAA 
0.0049436394 0.0000000000 0.0005768219 0.0024290355 0.0019873742 0.0000000000 
       TCTAC        TCTAG        TCTAT        TCTAN        TCTCA        TCTCC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TCTCG        TCTCT        TCTCN        TCTGA        TCTGC        TCTGG 
0.0000000000 0.0000000000 0.0000000000 0.0021556238 0.0021556238 0.0021556238 
       TCTGT        TCTGN        TCTTA        TCTTC        TCTTG        TCTTT 
0.0021556238 0.0021556238 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TCTTN        TCTNA        TCTNC        TCTNG        TCTNT        TCTNN 
0.0000000000 0.0005389060 0.0005389060 0.0005389060 0.0005389060 0.0005389060 
       TCNAA        TCNAC        TCNAG        TCNAT        TCNAN        TCNCA 
          NA           NA           NA           NA           NA           NA 
       TCNCC        TCNCG        TCNCT        TCNCN        TCNGA        TCNGC 
          NA           NA           NA           NA           NA           NA 
       TCNGG        TCNGT        TCNGN        TCNTA        TCNTC        TCNTG 
          NA           NA           NA           NA           NA           NA 
       TCNTT        TCNTN        TCNNA        TCNNC        TCNNG        TCNNT 
          NA           NA           NA           NA           NA           NA 
       TCNNN        TGAAA        TGAAC        TGAAG        TGAAT        TGAAN 
          NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TGACA        TGACC        TGACG        TGACT        TGACN        TGAGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0006253815 
       TGAGC        TGAGG        TGAGT        TGAGN        TGATA        TGATC 
0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0000000000 0.0000000000 
       TGATG        TGATT        TGATN        TGANA        TGANC        TGANG 
0.0000000000 0.0000000000 0.0000000000 0.0001563454 0.0001563454 0.0001563454 
       TGANT        TGANN        TGCAA        TGCAC        TGCAG        TGCAT 
0.0001563454 0.0001563454 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TGCAN        TGCCA        TGCCC        TGCCG        TGCCT        TGCCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TGCGA        TGCGC        TGCGG        TGCGT        TGCGN        TGCTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TGCTC        TGCTG        TGCTT        TGCTN        TGCNA        TGCNC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TGCNG        TGCNT        TGCNN        TGGAA        TGGAC        TGGAG 
0.0000000000 0.0000000000 0.0000000000 0.0018693956 0.0000000000 0.0023072876 
       TGGAT        TGGAN        TGGCA        TGGCC        TGGCG        TGGCT 
0.0000000000 0.0010441708 0.0025595886 0.0000000000 0.0000000000 0.0083306927 
       TGGCN        TGGGA        TGGGC        TGGGG        TGGGT        TGGGN 
0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 
       TGGTA        TGGTC        TGGTG        TGGTT        TGGTN        TGGNA 
0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0049436394 
       TGGNC        TGGNG        TGGNT        TGGNN        TGTAA        TGTAC 
0.0000000000 0.0005768219 0.0024290355 0.0019873742 0.0000000000 0.0000000000 
       TGTAG        TGTAT        TGTAN        TGTCA        TGTCC        TGTCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TGTCT        TGTCN        TGTGA        TGTGC        TGTGG        TGTGT 
0.0000000000 0.0000000000 0.0016522853 0.0016522853 0.0016522853 0.0016522853 
       TGTGN        TGTTA        TGTTC        TGTTG        TGTTT        TGTTN 
0.0016522853 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TGTNA        TGTNC        TGTNG        TGTNT        TGTNN        TGNAA 
0.0004130713 0.0004130713 0.0004130713 0.0004130713 0.0004130713           NA 
       TGNAC        TGNAG        TGNAT        TGNAN        TGNCA        TGNCC 
          NA           NA           NA           NA           NA           NA 
       TGNCG        TGNCT        TGNCN        TGNGA        TGNGC        TGNGG 
          NA           NA           NA           NA           NA           NA 
       TGNGT        TGNGN        TGNTA        TGNTC        TGNTG        TGNTT 
          NA           NA           NA           NA           NA           NA 
       TGNTN        TGNNA        TGNNC        TGNNG        TGNNT        TGNNN 
          NA           NA           NA           NA           NA           NA 
       TTAAA        TTAAC        TTAAG        TTAAT        TTAAN        TTACA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTACC        TTACG        TTACT        TTACN        TTAGA        TTAGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTAGG        TTAGT        TTAGN        TTATA        TTATC        TTATG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTATT        TTATN        TTANA        TTANC        TTANG        TTANT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTANN        TTCAA        TTCAC        TTCAG        TTCAT        TTCAN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTCCA        TTCCC        TTCCG        TTCCT        TTCCN        TTCGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTCGC        TTCGG        TTCGT        TTCGN        TTCTA        TTCTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTCTG        TTCTT        TTCTN        TTCNA        TTCNC        TTCNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTCNT        TTCNN        TTGAA        TTGAC        TTGAG        TTGAT 
0.0000000000 0.0000000000 0.0027781351 0.0000000000 0.0023072876 0.0000000000 
       TTGAN        TTGCA        TTGCC        TTGCG        TTGCT        TTGCN 
0.0012713557 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 
       TTGGA        TTGGC        TTGGG        TTGGT        TTGGN        TTGTA 
0.0005415140 0.0000000000 0.0000000000 0.0024857935 0.0007568269 0.0148040595 
       TTGTC        TTGTG        TTGTT        TTGTN        TTGNA        TTGNC 
0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0051708243 0.0000000000 
       TTGNG        TTGNT        TTGNN        TTTAA        TTTAC        TTTAG 
0.0005768219 0.0027041215 0.0021129419 0.0000000000 0.0000000000 0.0000000000 
       TTTAT        TTTAN        TTTCA        TTTCC        TTTCG        TTTCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTTCN        TTTGA        TTTGC        TTTGG        TTTGT        TTTGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTTTA        TTTTC        TTTTG        TTTTT        TTTTN        TTTNA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TTTNC        TTTNG        TTTNT        TTTNN        TTNAA        TTNAC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000           NA           NA 
       TTNAG        TTNAT        TTNAN        TTNCA        TTNCC        TTNCG 
          NA           NA           NA           NA           NA           NA 
       TTNCT        TTNCN        TTNGA        TTNGC        TTNGG        TTNGT 
          NA           NA           NA           NA           NA           NA 
       TTNGN        TTNTA        TTNTC        TTNTG        TTNTT        TTNTN 
          NA           NA           NA           NA           NA           NA 
       TTNNA        TTNNC        TTNNG        TTNNT        TTNNN        TNAAA 
          NA           NA           NA           NA           NA 0.0000000000 
       TNAAC        TNAAG        TNAAT        TNAAN        TNACA        TNACC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TNACG        TNACT        TNACN        TNAGA        TNAGC        TNAGG 
0.0000000000 0.0000000000 0.0000000000 0.0010554064 0.0010554064 0.0010554064 
       TNAGT        TNAGN        TNATA        TNATC        TNATG        TNATT 
0.0010554064 0.0010554064 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TNATN        TNANA        TNANC        TNANG        TNANT        TNANN 
0.0000000000 0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0002638516 
       TNCAA        TNCAC        TNCAG        TNCAT        TNCAN        TNCCA 
0.0027454851 0.0027454851 0.0027454851 0.0027454851 0.0027454851 0.0027454851 
       TNCCC        TNCCG        TNCCT        TNCCN        TNCGA        TNCGC 
0.0027454851 0.0027454851 0.0027454851 0.0027454851 0.0027454851 0.0027454851 
       TNCGG        TNCGT        TNCGN        TNCTA        TNCTC        TNCTG 
0.0027454851 0.0027454851 0.0027454851 0.0027454851 0.0027454851 0.0027454851 
       TNCTT        TNCTN        TNCNA        TNCNC        TNCNG        TNCNT 
0.0027454851 0.0027454851 0.0027454851 0.0027454851 0.0027454851 0.0027454851 
       TNCNN        TNGAA        TNGAC        TNGAG        TNGAT        TNGAN 
0.0027454851 0.0018693956 0.0000000000 0.0023072876 0.0000000000 0.0010441708 
       TNGCA        TNGCC        TNGCG        TNGCT        TNGCN        TNGGA 
0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 
       TNGGC        TNGGG        TNGGT        TNGGN        TNGTA        TNGTC 
0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 
       TNGTG        TNGTT        TNGTN        TNGNA        TNGNC        TNGNG 
0.0000000000 0.0000000000 0.0037010149 0.0049436394 0.0000000000 0.0005768219 
       TNGNT        TNGNN        TNTAA        TNTAC        TNTAG        TNTAT 
0.0024290355 0.0019873742 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TNTAN        TNTCA        TNTCC        TNTCG        TNTCT        TNTCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       TNTGA        TNTGC        TNTGG        TNTGT        TNTGN        TNTTA 
0.0009519773 0.0009519773 0.0009519773 0.0009519773 0.0009519773 0.0001450970 
       TNTTC        TNTTG        TNTTT        TNTTN        TNTNA        TNTNC 
0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0002742686 0.0002742686 
       TNTNG        TNTNT        TNTNN        TNNAA        TNNAC        TNNAG 
0.0002742686 0.0002742686 0.0002742686           NA           NA           NA 
       TNNAT        TNNAN        TNNCA        TNNCC        TNNCG        TNNCT 
          NA           NA           NA           NA           NA           NA 
       TNNCN        TNNGA        TNNGC        TNNGG        TNNGT        TNNGN 
          NA           NA           NA           NA           NA           NA 
       TNNTA        TNNTC        TNNTG        TNNTT        TNNTN        TNNNA 
          NA           NA           NA           NA           NA           NA 
       TNNNC        TNNNG        TNNNT        TNNNN        NAAAA        NAAAC 
          NA           NA           NA           NA 0.0000000000 0.0000000000 
       NAAAG        NAAAT        NAAAN        NAACA        NAACC        NAACG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NAACT        NAACN        NAAGA        NAAGC        NAAGG        NAAGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NAAGN        NAATA        NAATC        NAATG        NAATT        NAATN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NAANA        NAANC        NAANG        NAANT        NAANN        NACAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0028782383 
       NACAC        NACAG        NACAT        NACAN        NACCA        NACCC 
0.0028782383 0.0028782383 0.0028782383 0.0028782383 0.0028782383 0.0028782383 
       NACCG        NACCT        NACCN        NACGA        NACGC        NACGG 
0.0028782383 0.0028782383 0.0028782383 0.0028782383 0.0028782383 0.0028782383 
       NACGT        NACGN        NACTA        NACTC        NACTG        NACTT 
0.0028782383 0.0028782383 0.0028782383 0.0028782383 0.0028782383 0.0028782383 
       NACTN        NACNA        NACNC        NACNG        NACNT        NACNN 
0.0028782383 0.0028782383 0.0028782383 0.0028782383 0.0028782383 0.0028782383 
       NAGAA        NAGAC        NAGAG        NAGAT        NAGAN        NAGCA 
0.0016422107 0.0000000000 0.0023072876 0.0000000000 0.0009873746 0.0025595886 
       NAGCC        NAGCG        NAGCT        NAGCN        NAGGA        NAGGC 
0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 
       NAGGG        NAGGT        NAGGN        NAGTA        NAGTC        NAGTG 
0.0000000000 0.0008352771 0.0003441978 0.0148040595 0.0000000000 0.0000000000 
       NAGTT        NAGTN        NAGNA        NAGNC        NAGNG        NAGNT 
0.0000000000 0.0037010149 0.0048868432 0.0000000000 0.0005768219 0.0022914924 
       NAGNN        NATAA        NATAC        NATAG        NATAT        NATAN 
0.0019387894 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NATCA        NATCC        NATCG        NATCT        NATCN        NATGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NATGC        NATGG        NATGT        NATGN        NATTA        NATTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005803880 0.0005803880 
       NATTG        NATTT        NATTN        NATNA        NATNC        NATNG 
0.0005803880 0.0005803880 0.0005803880 0.0001450970 0.0001450970 0.0001450970 
       NATNT        NATNN        NANAA        NANAC        NANAG        NANAT 
0.0001450970 0.0001450970           NA           NA           NA           NA 
       NANAN        NANCA        NANCC        NANCG        NANCT        NANCN 
          NA           NA           NA           NA           NA           NA 
       NANGA        NANGC        NANGG        NANGT        NANGN        NANTA 
          NA           NA           NA           NA           NA           NA 
       NANTC        NANTG        NANTT        NANTN        NANNA        NANNC 
          NA           NA           NA           NA           NA           NA 
       NANNG        NANNT        NANNN        NCAAA        NCAAC        NCAAG 
          NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 
       NCAAT        NCAAN        NCACA        NCACC        NCACG        NCACT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NCACN        NCAGA        NCAGC        NCAGG        NCAGT        NCAGN 
0.0000000000 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0035962442 
       NCATA        NCATC        NCATG        NCATT        NCATN        NCANA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0008990611 
       NCANC        NCANG        NCANT        NCANN        NCCAA        NCCAC 
0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0011860551 0.0011860551 
       NCCAG        NCCAT        NCCAN        NCCCA        NCCCC        NCCCG 
0.0011860551 0.0004128226 0.0009927470 0.0007994389 0.0007994389 0.0007994389 
       NCCCT        NCCCN        NCCGA        NCCGC        NCCGG        NCCGT 
0.0007994389 0.0007994389 0.0007994389 0.0007994389 0.0007994389 0.0007994389 
       NCCGN        NCCTA        NCCTC        NCCTG        NCCTT        NCCTN 
0.0007994389 0.0004128226 0.0004128226 0.0004128226 0.0011860551 0.0006061307 
       NCCNA        NCCNC        NCCNG        NCCNT        NCCNN        NCGAA 
0.0007994389 0.0007994389 0.0007994389 0.0007994389 0.0007994389 0.0018693956 
       NCGAC        NCGAG        NCGAT        NCGAN        NCGCA        NCGCC 
0.0000000000 0.0023072876 0.0000000000 0.0010441708 0.0025595886 0.0000000000 
       NCGCG        NCGCT        NCGCN        NCGGA        NCGGC        NCGGG 
0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 
       NCGGT        NCGGN        NCGTA        NCGTC        NCGTG        NCGTT 
0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 
       NCGTN        NCGNA        NCGNC        NCGNG        NCGNT        NCGNN 
0.0037010149 0.0049436394 0.0000000000 0.0005768219 0.0024290355 0.0019873742 
       NCTAA        NCTAC        NCTAG        NCTAT        NCTAN        NCTCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NCTCC        NCTCG        NCTCT        NCTCN        NCTGA        NCTGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0021556238 0.0021556238 
       NCTGG        NCTGT        NCTGN        NCTTA        NCTTC        NCTTG 
0.0021556238 0.0021556238 0.0021556238 0.0000000000 0.0000000000 0.0000000000 
       NCTTT        NCTTN        NCTNA        NCTNC        NCTNG        NCTNT 
0.0000000000 0.0000000000 0.0005389060 0.0005389060 0.0005389060 0.0005389060 
       NCTNN        NCNAA        NCNAC        NCNAG        NCNAT        NCNAN 
0.0005389060           NA           NA           NA           NA           NA 
       NCNCA        NCNCC        NCNCG        NCNCT        NCNCN        NCNGA 
          NA           NA           NA           NA           NA           NA 
       NCNGC        NCNGG        NCNGT        NCNGN        NCNTA        NCNTC 
          NA           NA           NA           NA           NA           NA 
       NCNTG        NCNTT        NCNTN        NCNNA        NCNNC        NCNNG 
          NA           NA           NA           NA           NA           NA 
       NCNNT        NCNNN        NGAAA        NGAAC        NGAAG        NGAAT 
          NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NGAAN        NGACA        NGACC        NGACG        NGACT        NGACN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NGAGA        NGAGC        NGAGG        NGAGT        NGAGN        NGATA 
0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0000000000 
       NGATC        NGATG        NGATT        NGATN        NGANA        NGANC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0001563454 0.0001563454 
       NGANG        NGANT        NGANN        NGCAA        NGCAC        NGCAG 
0.0001563454 0.0001563454 0.0001563454 0.0018453454 0.0018453454 0.0018453454 
       NGCAT        NGCAN        NGCCA        NGCCC        NGCCG        NGCCT 
0.0018453454 0.0018453454 0.0018453454 0.0018453454 0.0018453454 0.0018453454 
       NGCCN        NGCGA        NGCGC        NGCGG        NGCGT        NGCGN 
0.0018453454 0.0018453454 0.0018453454 0.0018453454 0.0018453454 0.0018453454 
       NGCTA        NGCTC        NGCTG        NGCTT        NGCTN        NGCNA 
0.0018453454 0.0018453454 0.0018453454 0.0018453454 0.0018453454 0.0018453454 
       NGCNC        NGCNG        NGCNT        NGCNN        NGGAA        NGGAC 
0.0018453454 0.0018453454 0.0018453454 0.0018453454 0.0018693956 0.0000000000 
       NGGAG        NGGAT        NGGAN        NGGCA        NGGCC        NGGCG 
0.0023072876 0.0000000000 0.0010441708 0.0025595886 0.0000000000 0.0000000000 
       NGGCT        NGGCN        NGGGA        NGGGC        NGGGG        NGGGT 
0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 
       NGGGN        NGGTA        NGGTC        NGGTG        NGGTT        NGGTN 
0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 
       NGGNA        NGGNC        NGGNG        NGGNT        NGGNN        NGTAA 
0.0049436394 0.0000000000 0.0005768219 0.0024290355 0.0019873742 0.0000000000 
       NGTAC        NGTAG        NGTAT        NGTAN        NGTCA        NGTCC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NGTCG        NGTCT        NGTCN        NGTGA        NGTGC        NGTGG 
0.0000000000 0.0000000000 0.0000000000 0.0016522853 0.0016522853 0.0016522853 
       NGTGT        NGTGN        NGTTA        NGTTC        NGTTG        NGTTT 
0.0016522853 0.0016522853 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NGTTN        NGTNA        NGTNC        NGTNG        NGTNT        NGTNN 
0.0000000000 0.0004130713 0.0004130713 0.0004130713 0.0004130713 0.0004130713 
       NGNAA        NGNAC        NGNAG        NGNAT        NGNAN        NGNCA 
          NA           NA           NA           NA           NA           NA 
       NGNCC        NGNCG        NGNCT        NGNCN        NGNGA        NGNGC 
          NA           NA           NA           NA           NA           NA 
       NGNGG        NGNGT        NGNGN        NGNTA        NGNTC        NGNTG 
          NA           NA           NA           NA           NA           NA 
       NGNTT        NGNTN        NGNNA        NGNNC        NGNNG        NGNNT 
          NA           NA           NA           NA           NA           NA 
       NGNNN        NTAAA        NTAAC        NTAAG        NTAAT        NTAAN 
          NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTACA        NTACC        NTACG        NTACT        NTACN        NTAGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTAGC        NTAGG        NTAGT        NTAGN        NTATA        NTATC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTATG        NTATT        NTATN        NTANA        NTANC        NTANG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTANT        NTANN        NTCAA        NTCAC        NTCAG        NTCAT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTCAN        NTCCA        NTCCC        NTCCG        NTCCT        NTCCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTCGA        NTCGC        NTCGG        NTCGT        NTCGN        NTCTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTCTC        NTCTG        NTCTT        NTCTN        NTCNA        NTCNC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTCNG        NTCNT        NTCNN        NTGAA        NTGAC        NTGAG 
0.0000000000 0.0000000000 0.0000000000 0.0020965805 0.0000000000 0.0023072876 
       NTGAT        NTGAN        NTGCA        NTGCC        NTGCG        NTGCT 
0.0000000000 0.0011009670 0.0025595886 0.0000000000 0.0000000000 0.0083306927 
       NTGCN        NTGGA        NTGGC        NTGGG        NTGGT        NTGGN 
0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0019356213 0.0006192838 
       NTGTA        NTGTC        NTGTG        NTGTT        NTGTN        NTGNA 
0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0050004356 
       NTGNC        NTGNG        NTGNT        NTGNN        NTTAA        NTTAC 
0.0000000000 0.0005768219 0.0025665785 0.0020359590 0.0000000000 0.0000000000 
       NTTAG        NTTAT        NTTAN        NTTCA        NTTCC        NTTCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTTCT        NTTCN        NTTGA        NTTGC        NTTGG        NTTGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTTGN        NTTTA        NTTTC        NTTTG        NTTTT        NTTTN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NTTNA        NTTNC        NTTNG        NTTNT        NTTNN        NTNAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000           NA 
       NTNAC        NTNAG        NTNAT        NTNAN        NTNCA        NTNCC 
          NA           NA           NA           NA           NA           NA 
       NTNCG        NTNCT        NTNCN        NTNGA        NTNGC        NTNGG 
          NA           NA           NA           NA           NA           NA 
       NTNGT        NTNGN        NTNTA        NTNTC        NTNTG        NTNTT 
          NA           NA           NA           NA           NA           NA 
       NTNTN        NTNNA        NTNNC        NTNNG        NTNNT        NTNNN 
          NA           NA           NA           NA           NA           NA 
       NNAAA        NNAAC        NNAAG        NNAAT        NNAAN        NNACA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NNACC        NNACG        NNACT        NNACN        NNAGA        NNAGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0010554064 0.0010554064 
       NNAGG        NNAGT        NNAGN        NNATA        NNATC        NNATG 
0.0010554064 0.0010554064 0.0010554064 0.0000000000 0.0000000000 0.0000000000 
       NNATT        NNATN        NNANA        NNANC        NNANG        NNANT 
0.0000000000 0.0000000000 0.0002638516 0.0002638516 0.0002638516 0.0002638516 
       NNANN        NNCAA        NNCAC        NNCAG        NNCAT        NNCAN 
0.0002638516 0.0014774097 0.0014774097 0.0014774097 0.0012841016 0.0014290827 
       NNCCA        NNCCC        NNCCG        NNCCT        NNCCN        NNCGA 
0.0013807556 0.0013807556 0.0013807556 0.0013807556 0.0013807556 0.0013807556 
       NNCGC        NNCGG        NNCGT        NNCGN        NNCTA        NNCTC 
0.0013807556 0.0013807556 0.0013807556 0.0013807556 0.0012841016 0.0012841016 
       NNCTG        NNCTT        NNCTN        NNCNA        NNCNC        NNCNG 
0.0012841016 0.0014774097 0.0013324286 0.0013807556 0.0013807556 0.0013807556 
       NNCNT        NNCNN        NNGAA        NNGAC        NNGAG        NNGAT 
0.0013807556 0.0013807556 0.0018693956 0.0000000000 0.0023072876 0.0000000000 
       NNGAN        NNGCA        NNGCC        NNGCG        NNGCT        NNGCN 
0.0010441708 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 
       NNGGA        NNGGC        NNGGG        NNGGT        NNGGN        NNGTA 
0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 
       NNGTC        NNGTG        NNGTT        NNGTN        NNGNA        NNGNC 
0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0049436394 0.0000000000 
       NNGNG        NNGNT        NNGNN        NNTAA        NNTAC        NNTAG 
0.0005768219 0.0024290355 0.0019873742 0.0000000000 0.0000000000 0.0000000000 
       NNTAT        NNTAN        NNTCA        NNTCC        NNTCG        NNTCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       NNTCN        NNTGA        NNTGC        NNTGG        NNTGT        NNTGN 
0.0000000000 0.0009519773 0.0009519773 0.0009519773 0.0009519773 0.0009519773 
       NNTTA        NNTTC        NNTTG        NNTTT        NNTTN        NNTNA 
0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0002742686 
       NNTNC        NNTNG        NNTNT        NNTNN        NNNAA        NNNAC 
0.0002742686 0.0002742686 0.0002742686 0.0002742686           NA           NA 
       NNNAG        NNNAT        NNNAN        NNNCA        NNNCC        NNNCG 
          NA           NA           NA           NA           NA           NA 
       NNNCT        NNNCN        NNNGA        NNNGC        NNNGG        NNNGT 
          NA           NA           NA           NA           NA           NA 
       NNNGN        NNNTA        NNNTC        NNNTG        NNNTT        NNNTN 
          NA           NA           NA           NA           NA           NA 
       NNNNA        NNNNC        NNNNG        NNNNT        NNNNN 
          NA           NA           NA           NA           NA 

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






