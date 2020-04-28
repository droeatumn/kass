Submit-block ::= {
  contact {
    contact {
      name name {
        last "Clark",
        first "Karen"
      },
      affil std {
        affil "NIH",
        div "NCBI",
        city "Bethesda",
        sub "MD",
        country "USA",
        street "Center Dr.",
        email "kclark@ncbi.nlm.nih.gov",
        fax "301-444-5555",
        phone "301-435-5936",
        postal-code "20872"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "Clark",
            first "Karen",
            initials "K.L.",
            suffix ""
          }
        },
        {
          name name {
            last "Smith",
            first "Thomas",
            initials "T.",
            suffix ""
          }
        }
      },
      affil std {
        affil "NIH",
        div "NCBI",
        city "Bethesda",
        sub "MD",
        country "USA",
        street "Center Dr.",
        postal-code "20872"
      }
    }
  },
  subtype new
}

Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "Clark",
              first "Karen",
              initials "K.L.",
              suffix ""
            }
          },
          {
            name name {
              last "Smith",
              first "Thomas",
              initials "T.",
              suffix ""
            }
          },
          {
            name name {
              last "Jones",
              first "Nicolas",
              initials "N.F.",
              suffix ""
            }
          },
          {
            name name {
              last "Wesson",
              first "Huntington",
              initials "H.F.",
              suffix ""
            }
          }
        }
      },
      title "Draft sequence of the ABC genome"
    }
  }
}

Seqdesc ::= user {
  type str "DBLink",
  data {
    {
      label str "BioProject",
      num 1,
      data strs {
        "PRJNA33011"
      }
    }
  }
}

