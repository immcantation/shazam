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
       AAAAA        AAAAC        AAAAG        AAAAT        AAAAN        AAACA        AAACC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAACG        AAACT        AAACN        AAAGA        AAAGC        AAAGG        AAAGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAAGN        AAATA        AAATC        AAATG        AAATT        AAATN        AAANA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AAANC        AAANG        AAANT        AAANN        AACAA        AACAC        AACAG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACAT        AACAN        AACCA        AACCC        AACCG        AACCT        AACCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACGA        AACGC        AACGG        AACGT        AACGN        AACTA        AACTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACTG        AACTT        AACTN        AACNA        AACNC        AACNG        AACNT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AACNN        AAGAA        AAGAC        AAGAG        AAGAT        AAGAN        AAGCA 
0.0000000000 0.0027781351 0.0000000000 0.0023072876 0.0000000000 0.0012713557 0.0025595886 
       AAGCC        AAGCG        AAGCT        AAGCN        AAGGA        AAGGC        AAGGG 
0.0000000000 0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 
       AAGGT        AAGGN        AAGTA        AAGTC        AAGTG        AAGTT        AAGTN 
0.0002851050 0.0002066547 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 
       AAGNA        AAGNC        AAGNG        AAGNT        AAGNN        AATAA        AATAC 
0.0051708243 0.0000000000 0.0005768219 0.0021539494 0.0019753989 0.0000000000 0.0000000000 
       AATAG        AATAT        AATAN        AATCA        AATCC        AATCG        AATCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AATCN        AATGA        AATGC        AATGG        AATGT        AATGN        AATTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005803880 
       AATTC        AATTG        AATTT        AATTN        AATNA        AATNC        AATNG 
0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0001450970 0.0001450970 0.0001450970 
       AATNT        AATNN        AANAA        AANAC        AANAG        AANAT        AANAN 
0.0001450970 0.0001450970           NA           NA           NA           NA           NA 
       AANCA        AANCC        AANCG        AANCT        AANCN        AANGA        AANGC 
          NA           NA           NA           NA           NA           NA           NA 
       AANGG        AANGT        AANGN        AANTA        AANTC        AANTG        AANTT 
          NA           NA           NA           NA           NA           NA           NA 
       AANTN        AANNA        AANNC        AANNG        AANNT        AANNN        ACAAA 
          NA           NA           NA           NA           NA           NA 0.0000000000 
       ACAAC        ACAAG        ACAAT        ACAAN        ACACA        ACACC        ACACG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACACT        ACACN        ACAGA        ACAGC        ACAGG        ACAGT        ACAGN 
0.0000000000 0.0000000000 0.0035962442 0.0035962442 0.0035962442 0.0035962442 0.0035962442 
       ACATA        ACATC        ACATG        ACATT        ACATN        ACANA        ACANC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0008990611 0.0008990611 
       ACANG        ACANT        ACANN        ACCAA        ACCAC        ACCAG        ACCAT 
0.0008990611 0.0008990611 0.0008990611 0.0043115900 0.0043115900 0.0043115900 0.0012186597 
       ACCAN        ACCCA        ACCCC        ACCCG        ACCCT        ACCCN        ACCGA 
0.0035383574 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 
       ACCGC        ACCGG        ACCGT        ACCGN        ACCTA        ACCTC        ACCTG 
0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0012186597 0.0012186597 0.0012186597 
       ACCTT        ACCTN        ACCNA        ACCNC        ACCNG        ACCNT        ACCNN 
0.0043115900 0.0019918923 0.0027651248 0.0027651248 0.0027651248 0.0027651248 0.0027651248 
       ACGAA        ACGAC        ACGAG        ACGAT        ACGAN        ACGCA        ACGCC 
0.0009606561 0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0025595886 0.0000000000 
       ACGCG        ACGCT        ACGCN        ACGGA        ACGGC        ACGGG        ACGGT 
0.0000000000 0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 
       ACGGN        ACGTA        ACGTC        ACGTG        ACGTT        ACGTN        ACGNA 
0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0047164545 
       ACGNC        ACGNG        ACGNT        ACGNN        ACTAA        ACTAC        ACTAG 
0.0000000000 0.0005768219 0.0024290355 0.0019305780 0.0000000000 0.0000000000 0.0000000000 
       ACTAT        ACTAN        ACTCA        ACTCC        ACTCG        ACTCT        ACTCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ACTGA        ACTGC        ACTGG        ACTGT        ACTGN        ACTTA        ACTTC 
0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0000000000 0.0000000000 
       ACTTG        ACTTT        ACTTN        ACTNA        ACTNC        ACTNG        ACTNT 
0.0000000000 0.0000000000 0.0000000000 0.0005389060 0.0005389060 0.0005389060 0.0005389060 
       ACTNN        ACNAA        ACNAC        ACNAG        ACNAT        ACNAN        ACNCA 
0.0005389060           NA           NA           NA           NA           NA           NA 
       ACNCC        ACNCG        ACNCT        ACNCN        ACNGA        ACNGC        ACNGG 
          NA           NA           NA           NA           NA           NA           NA 
       ACNGT        ACNGN        ACNTA        ACNTC        ACNTG        ACNTT        ACNTN 
          NA           NA           NA           NA           NA           NA           NA 
       ACNNA        ACNNC        ACNNG        ACNNT        ACNNN        AGAAA        AGAAC 
          NA           NA           NA           NA           NA 0.0000000000 0.0000000000 
       AGAAG        AGAAT        AGAAN        AGACA        AGACC        AGACG        AGACT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AGACN        AGAGA        AGAGC        AGAGG        AGAGT        AGAGN        AGATA 
0.0000000000 0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0006253815 0.0000000000 
       AGATC        AGATG        AGATT        AGATN        AGANA        AGANC        AGANG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0001563454 0.0001563454 0.0001563454 
       AGANT        AGANN        AGCAA        AGCAC        AGCAG        AGCAT        AGCAN 
0.0001563454 0.0001563454 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCCA        AGCCC        AGCCG        AGCCT        AGCCN        AGCGA        AGCGC 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCGG        AGCGT        AGCGN        AGCTA        AGCTC        AGCTG        AGCTT 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 
       AGCTN        AGCNA        AGCNC        AGCNG        AGCNT        AGCNN        AGGAA 
0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0073813816 0.0009606561 
       AGGAC        AGGAG        AGGAT        AGGAN        AGGCA        AGGCC        AGGCG 
0.0000000000 0.0023072876 0.0000000000 0.0008169859 0.0025595886 0.0000000000 0.0000000000 
       AGGCT        AGGCN        AGGGA        AGGGC        AGGGG        AGGGT        AGGGN 
0.0083306927 0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 
       AGGTA        AGGTC        AGGTG        AGGTT        AGGTN        AGGNA        AGGNC 
0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0047164545 0.0000000000 
       AGGNG        AGGNT        AGGNN        AGTAA        AGTAC        AGTAG        AGTAT 
0.0005768219 0.0024290355 0.0019305780 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       AGTAN        AGTCA        AGTCC        AGTCG        AGTCT        AGTCN        AGTGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0016522853 
       AGTGC        AGTGG        AGTGT        AGTGN        AGTTA        AGTTC        AGTTG 
0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0000000000 0.0000000000 0.0000000000 
       AGTTT        AGTTN        AGTNA        AGTNC        AGTNG        AGTNT        AGTNN 
0.0000000000 0.0000000000 0.0004130713 0.0004130713 0.0004130713 0.0004130713 0.0004130713 
       AGNAA        AGNAC        AGNAG        AGNAT        AGNAN        AGNCA        AGNCC 
          NA           NA           NA           NA           NA           NA           NA 
       AGNCG        AGNCT        AGNCN        AGNGA        AGNGC        AGNGG        AGNGT 
          NA           NA           NA           NA           NA           NA           NA 
       AGNGN        AGNTA        AGNTC        AGNTG        AGNTT        AGNTN        AGNNA 
          NA           NA           NA           NA           NA           NA           NA 
       AGNNC        AGNNG        AGNNT        AGNNN        ATAAA        ATAAC        ATAAG 
          NA           NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 
       ATAAT        ATAAN        ATACA        ATACC        ATACG        ATACT        ATACN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATAGA        ATAGC        ATAGG        ATAGT        ATAGN        ATATA        ATATC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATATG        ATATT        ATATN        ATANA        ATANC        ATANG        ATANT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATANN        ATCAA        ATCAC        ATCAG        ATCAT        ATCAN        ATCCA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCCC        ATCCG        ATCCT        ATCCN        ATCGA        ATCGC        ATCGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCGT        ATCGN        ATCTA        ATCTC        ATCTG        ATCTT        ATCTN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATCNA        ATCNC        ATCNG        ATCNT        ATCNN        ATGAA        ATGAC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0018693956 0.0000000000 
       ATGAG        ATGAT        ATGAN        ATGCA        ATGCC        ATGCG        ATGCT 
0.0023072876 0.0000000000 0.0010441708 0.0025595886 0.0000000000 0.0000000000 0.0083306927 
       ATGCN        ATGGA        ATGGC        ATGGG        ATGGT        ATGGN        ATGTA 
0.0027225703 0.0005415140 0.0000000000 0.0000000000 0.0024857935 0.0007568269 0.0148040595 
       ATGTC        ATGTG        ATGTT        ATGTN        ATGNA        ATGNC        ATGNG 
0.0000000000 0.0000000000 0.0000000000 0.0037010149 0.0049436394 0.0000000000 0.0005768219 
       ATGNT        ATGNN        ATTAA        ATTAC        ATTAG        ATTAT        ATTAN 
0.0027041215 0.0020561457 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTCA        ATTCC        ATTCG        ATTCT        ATTCN        ATTGA        ATTGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTGG        ATTGT        ATTGN        ATTTA        ATTTC        ATTTG        ATTTT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ATTTN        ATTNA        ATTNC        ATTNG        ATTNT        ATTNN        ATNAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000           NA 
       ATNAC        ATNAG        ATNAT        ATNAN        ATNCA        ATNCC        ATNCG 
          NA           NA           NA           NA           NA           NA           NA 
       ATNCT        ATNCN        ATNGA        ATNGC        ATNGG        ATNGT        ATNGN 
          NA           NA           NA           NA           NA           NA           NA 
       ATNTA        ATNTC        ATNTG        ATNTT        ATNTN        ATNNA        ATNNC 
          NA           NA           NA           NA           NA           NA           NA 
       ATNNG        ATNNT        ATNNN        ANAAA        ANAAC        ANAAG        ANAAT 
          NA           NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANAAN        ANACA        ANACC        ANACG        ANACT        ANACN        ANAGA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0010554064 
       ANAGC        ANAGG        ANAGT        ANAGN        ANATA        ANATC        ANATG 
0.0010554064 0.0010554064 0.0010554064 0.0010554064 0.0000000000 0.0000000000 0.0000000000 
       ANATT        ANATN        ANANA        ANANC        ANANG        ANANT        ANANN 
0.0000000000 0.0000000000 0.0002638516 0.0002638516 0.0002638516 0.0002638516 0.0002638516 
       ANCAA        ANCAC        ANCAG        ANCAT        ANCAN        ANCCA        ANCCC 
0.0029232429 0.0029232429 0.0029232429 0.0021500103 0.0027299347 0.0025366266 0.0025366266 
       ANCCG        ANCCT        ANCCN        ANCGA        ANCGC        ANCGG        ANCGT 
0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0025366266 
       ANCGN        ANCTA        ANCTC        ANCTG        ANCTT        ANCTN        ANCNA 
0.0025366266 0.0021500103 0.0021500103 0.0021500103 0.0029232429 0.0023433185 0.0025366266 
       ANCNC        ANCNG        ANCNT        ANCNN        ANGAA        ANGAC        ANGAG 
0.0025366266 0.0025366266 0.0025366266 0.0025366266 0.0016422107 0.0000000000 0.0023072876 
       ANGAT        ANGAN        ANGCA        ANGCC        ANGCG        ANGCT        ANGCN 
0.0000000000 0.0009873746 0.0025595886 0.0000000000 0.0000000000 0.0083306927 0.0027225703 
       ANGGA        ANGGC        ANGGG        ANGGT        ANGGN        ANGTA        ANGTC 
0.0005415140 0.0000000000 0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 
       ANGTG        ANGTT        ANGTN        ANGNA        ANGNC        ANGNG        ANGNT 
0.0000000000 0.0000000000 0.0037010149 0.0048868432 0.0000000000 0.0005768219 0.0024290355 
       ANGNN        ANTAA        ANTAC        ANTAG        ANTAT        ANTAN        ANTCA 
0.0019731751 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       ANTCC        ANTCG        ANTCT        ANTCN        ANTGA        ANTGC        ANTGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0009519773 0.0009519773 0.0009519773 
       ANTGT        ANTGN        ANTTA        ANTTC        ANTTG        ANTTT        ANTTN 
0.0009519773 0.0009519773 0.0001450970 0.0001450970 0.0001450970 0.0001450970 0.0001450970 
       ANTNA        ANTNC        ANTNG        ANTNT        ANTNN        ANNAA        ANNAC 
0.0002742686 0.0002742686 0.0002742686 0.0002742686 0.0002742686           NA           NA 
       ANNAG        ANNAT        ANNAN        ANNCA        ANNCC        ANNCG        ANNCT 
          NA           NA           NA           NA           NA           NA           NA 
       ANNCN        ANNGA        ANNGC        ANNGG        ANNGT        ANNGN        ANNTA 
          NA           NA           NA           NA           NA           NA           NA 
       ANNTC        ANNTG        ANNTT        ANNTN        ANNNA        ANNNC        ANNNG 
          NA           NA           NA           NA           NA           NA           NA 
       ANNNT        ANNNN        CAAAA        CAAAC        CAAAG        CAAAT        CAAAN 
          NA           NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAACA        CAACC        CAACG        CAACT        CAACN        CAAGA        CAAGC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAAGG        CAAGT        CAAGN        CAATA        CAATC        CAATG        CAATT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CAATN        CAANA        CAANC        CAANG        CAANT        CAANN        CACAA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACAC        CACAG        CACAT        CACAN        CACCA        CACCC        CACCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACCT        CACCN        CACGA        CACGC        CACGG        CACGT        CACGN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACTA        CACTC        CACTG        CACTT        CACTN        CACNA        CACNC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CACNG        CACNT        CACNN        CAGAA        CAGAC        CAGAG        CAGAT 
0.0000000000 0.0000000000 0.0000000000 0.0009606561 0.0000000000 0.0023072876 0.0000000000 
       CAGAN        CAGCA        CAGCC        CAGCG        CAGCT        CAGCN        CAGGA 
0.0008169859 0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 
       CAGGC        CAGGG        CAGGT        CAGGN        CAGTA        CAGTC        CAGTG 
0.0000000000 0.0000000000 0.0024857935 0.0007568269 0.0148040595 0.0000000000 0.0000000000 
       CAGTT        CAGTN        CAGNA        CAGNC        CAGNG        CAGNT        CAGNN 
0.0000000000 0.0037010149 0.0046447092 0.0000000000 0.0005768219 0.0027041215 0.0019814132 
       CATAA        CATAC        CATAG        CATAT        CATAN        CATCA        CATCC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CATCG        CATCT        CATCN        CATGA        CATGC        CATGG        CATGT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CATGN        CATTA        CATTC        CATTG        CATTT        CATTN        CATNA 
0.0000000000 0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0005803880 0.0001450970 
       CATNC        CATNG        CATNT        CATNN        CANAA        CANAC        CANAG 
0.0001450970 0.0001450970 0.0001450970 0.0001450970           NA           NA           NA 
       CANAT        CANAN        CANCA        CANCC        CANCG        CANCT        CANCN 
          NA           NA           NA           NA           NA           NA           NA 
       CANGA        CANGC        CANGG        CANGT        CANGN        CANTA        CANTC 
          NA           NA           NA           NA           NA           NA           NA 
       CANTG        CANTT        CANTN        CANNA        CANNC        CANNG        CANNT 
          NA           NA           NA           NA           NA           NA           NA 
       CANNN        CCAAA        CCAAC        CCAAG        CCAAT        CCAAN        CCACA 
          NA 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCACC        CCACG        CCACT        CCACN        CCAGA        CCAGC        CCAGG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0035962442 0.0035962442 0.0035962442 
       CCAGT        CCAGN        CCATA        CCATC        CCATG        CCATT        CCATN 
0.0035962442 0.0035962442 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCANA        CCANC        CCANG        CCANT        CCANN        CCCAA        CCCAC 
0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0008990611 0.0000000000 0.0000000000 
       CCCAG        CCCAT        CCCAN        CCCCA        CCCCC        CCCCG        CCCCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCCN        CCCGA        CCCGC        CCCGG        CCCGT        CCCGN        CCCTA 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCTC        CCCTG        CCCTT        CCCTN        CCCNA        CCCNC        CCCNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCCNT        CCCNN        CCGAA        CCGAC        CCGAG        CCGAT        CCGAN 
0.0000000000 0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 0.0010441708 
       CCGCA        CCGCC        CCGCG        CCGCT        CCGCN        CCGGA        CCGGC 
0.0022726072 0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 0.0000000000 
       CCGGG        CCGGT        CCGGN        CCGTA        CCGTC        CCGTG        CCGTT 
0.0000000000 0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 
       CCGTN        CCGNA        CCGNC        CCGNG        CCGNT        CCGNN        CCTAA 
0.0037010149 0.0048718941 0.0000000000 0.0005768219 0.0024290355 0.0019694379 0.0000000000 
       CCTAC        CCTAG        CCTAT        CCTAN        CCTCA        CCTCC        CCTCG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CCTCT        CCTCN        CCTGA        CCTGC        CCTGG        CCTGT        CCTGN 
0.0000000000 0.0000000000 0.0021556238 0.0021556238 0.0021556238 0.0021556238 0.0021556238 
       CCTTA        CCTTC        CCTTG        CCTTT        CCTTN        CCTNA        CCTNC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005389060 0.0005389060 
       CCTNG        CCTNT        CCTNN        CCNAA        CCNAC        CCNAG        CCNAT 
0.0005389060 0.0005389060 0.0005389060           NA           NA           NA           NA 
       CCNAN        CCNCA        CCNCC        CCNCG        CCNCT        CCNCN        CCNGA 
          NA           NA           NA           NA           NA           NA           NA 
       CCNGC        CCNGG        CCNGT        CCNGN        CCNTA        CCNTC        CCNTG 
          NA           NA           NA           NA           NA           NA           NA 
       CCNTT        CCNTN        CCNNA        CCNNC        CCNNG        CCNNT        CCNNN 
          NA           NA           NA           NA           NA           NA           NA 
       CGAAA        CGAAC        CGAAG        CGAAT        CGAAN        CGACA        CGACC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGACG        CGACT        CGACN        CGAGA        CGAGC        CGAGG        CGAGT 
0.0000000000 0.0000000000 0.0000000000 0.0006253815 0.0006253815 0.0006253815 0.0006253815 
       CGAGN        CGATA        CGATC        CGATG        CGATT        CGATN        CGANA 
0.0006253815 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0001563454 
       CGANC        CGANG        CGANT        CGANN        CGCAA        CGCAC        CGCAG 
0.0001563454 0.0001563454 0.0001563454 0.0001563454 0.0000000000 0.0000000000 0.0000000000 
       CGCAT        CGCAN        CGCCA        CGCCC        CGCCG        CGCCT        CGCCN 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCGA        CGCGC        CGCGG        CGCGT        CGCGN        CGCTA        CGCTC 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCTG        CGCTT        CGCTN        CGCNA        CGCNC        CGCNG        CGCNT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGCNN        CGGAA        CGGAC        CGGAG        CGGAT        CGGAN        CGGCA 
0.0000000000 0.0018693956 0.0000000000 0.0023072876 0.0000000000 0.0010441708 0.0022726072 
       CGGCC        CGGCG        CGGCT        CGGCN        CGGGA        CGGGC        CGGGG 
0.0000000000 0.0000000000 0.0083306927 0.0026508250 0.0005415140 0.0000000000 0.0000000000 
       CGGGT        CGGGN        CGGTA        CGGTC        CGGTG        CGGTT        CGGTN 
0.0013854492 0.0004817408 0.0148040595 0.0000000000 0.0000000000 0.0000000000 0.0037010149 
       CGGNA        CGGNC        CGGNG        CGGNT        CGGNN        CGTAA        CGTAC 
0.0048718941 0.0000000000 0.0005768219 0.0024290355 0.0019694379 0.0000000000 0.0000000000 
       CGTAG        CGTAT        CGTAN        CGTCA        CGTCC        CGTCG        CGTCT 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 
       CGTCN        CGTGA        CGTGC        CGTGG        CGTGT        CGTGN        CGTTA 
0.0000000000 0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0016522853 0.0000000000 
       CGTTC        CGTTG        CGTTT        CGTTN        CGTNA        CGTNC        CGTNG 
0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0004130713 0.0004130713 0.0004130713 
       CGTNT        CGTNN        CGNAA        CGNAC        CGNAG        CGNAT        CGNAN 
0.0004130713 0.0004130713           NA           NA           NA           NA           NA 
       CGNCA        CGNCC        CGNCG        CGNCT        CGNCN        CGNGA        CGNGC 
          NA           NA           NA           NA           NA           NA           NA 
       CGNGG        CGNGT        CGNGN        CGNTA        CGNTC        CGNTG        CGNTT 
          NA           NA           NA           NA           NA           NA           NA 
       CGNTN        CGNNA        CGNNC        CGNNG        CGNNT        CGNNN 
          NA           NA           NA           NA           NA           NA 
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






