sd_section("Overview",
           "Package overview",
           c("shazam"))

sd_section("Targeting models",
           "Targeting models",
           c("createTargetingModel", 
             "plotMutability",
             "createSubstitutionMatrix",
             "createMutabilityMatrix",
             "createTargetingMatrix",
             "extendSubstitutionMatrix",
             "extendMutabilityMatrix",
             "rescaleMutability",
             "TargetingModel-class",
             "M1NDistance",
             "HS1FDistance",
             "HS5FModel",
             "U5NModel"))

sd_section("Mutational profiling",
           "Mutational profiling",
           c("collapseByClone", 
             "calcDBObservedMutations", 
             "calcDBExpectedMutations",
             "calcObservedMutations",
             "calcExpectedMutations",
             "MutationDefinition-class",
             "RegionDefinition-class",
             "createMutationDefinition",
             "createRegionDefinition",
             "MUTATION_SCHEMES",
             "IMGT_SCHEMES"))

sd_section("Selection analysis",
           "Selection analysis",
           c("calcBaseline",
             "groupBaseline",
             "plotBaselineDensity",
             "plotBaselineSummary",
             "createBaseline",
             "getBaselineStats",
             "summarizeBaseline",
             "Baseline-class"))

sd_section("Distance profiling",
           "Distance profiling",
           c("distToNearest", 
             "calcTargetingDistance",
             "writeTargetingDistance"))

