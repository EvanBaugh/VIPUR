Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "psi",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "Query_1"
          },
          descr {
            title "2C35_chain_A",
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">2C35_chain_A"
                }
              }
            }
          },
          inst {
            repr raw,
            mol aa,
            length 129,
            seq-data ncbieaa "EEDASQLIFPKEFETAETLLNSEVHMLLEHRKQQNESAEDEQELSEVF
MKTLNYTARFSRFKNRETIASVRSLLLQKKLHKFELACLANLCPETAEESKALIPSLEGRFEDEELQQILDDIQTKRS
FQY"
          }
        }
      }
    },
    subject database "/home/evan/bio/databases/blast/nr/nr",
    algorithm-options {
      {
        name "InclusionThreshold",
        value real { 1, 10, -3 }
      },
      {
        name "PseudoCountWeight",
        value integer 2
      },
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 0 }
      },
      {
        name "MaskAtHash",
        value boolean FALSE
      },
      {
        name "SegFiltering",
        value boolean FALSE
      },
      {
        name "WordThreshold",
        value integer 11
      },
      {
        name "WindowSize",
        value integer 40
      },
      {
        name "HitlistSize",
        value integer 300
      },
      {
        name "CompositionBasedStats",
        value integer 1
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      },
      {
        name "GapTrigger",
        value real { 22, 10, 0 }
      }
    }
  }
}
