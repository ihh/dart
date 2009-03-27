#ifndef TRANSDUCER_KEYWORDS_INCLUDED
#define TRANSDUCER_KEYWORDS_INCLUDED


// keywords for transducer S-expressions

#define TSEXPR_TOKEN      "token"
#define TSEXPR_BRANCH     "branch"

#define TSEXPR_TRANSDUCER "transducer"
#define TSEXPR_COMPOSITE  "composite-transducer"
#define TSEXPR_ELIMINATED "acyclic-composite-transducer"

#define TSEXPR_VALUE      "value"
#define TSEXPR_BITVALUE   "bit-value"

#define TSEXPR_SEQUENCE   "sequence"
#define TSEXPR_PROFILE    "profile"
#define TSEXPR_BITPROFILE "bit-profile"

#define TSEXPR_RECONS     "reconstruct"

#define TSEXPR_FORMAT     "composite-name-format"
#define TSEXPR_PREFIX     "token-prefix"

#define TSEXPR_STATE      "state"
#define TSEXPR_NAME       "name"
#define TSEXPR_TYPE       "type"

#define TSEXPR_OLD_STATE  "old-state"
#define TSEXPR_BAND_COEFF "banding-coefficient"

#define TSEXPR_START      "start"
#define TSEXPR_END        "end"
#define TSEXPR_WAIT       "wait"
#define TSEXPR_MATCH      "match"
#define TSEXPR_DELETE     "delete"
#define TSEXPR_INSERT     "insert"
#define TSEXPR_EMIT       "emit"

#define TSEXPR_TRANSITION "transition"
#define TSEXPR_FROM       "from"
#define TSEXPR_TO         "to"
#define TSEXPR_LABEL      "label"
#define TSEXPR_SUM        "sum"

#define TSEXPR_PSCORE     "peeled-bits"
#define TSEXPR_FULL_SCORE "path-bits"
#define TSEXPR_OLD_SCORE  "old-path-bits"
#define TSEXPR_VSCORE     "viterbi-bits"
#define TSEXPR_FSCORE     "forward-bits"

#define TSEXPR_VPATH      "viterbi-path"
#define TSEXPR_FPATH      "forward-path"
#define TSEXPR_OPTACC     "optimal-accuracy-path"

#define TSEXPR_POSTEXPECT "expected-count"

#define TSEXPR_UNCONS     "unconstrained-subtree"
#define TSEXPR_CONS       "constrained-subtree"

#define TSEXPR_ALIGNMENT  "alignment"
#define TSEXPR_ALIGN_COL  "column"
#define TSEXPR_ID         "id"

#define TSEXPR_PEELED     "peeled-composition"

#define TSEXPR_SCHEDULE   "proposal-schedule"
#define TSEXPR_PROPOSAL   "proposal-composition"
#define TSEXPR_INVERSE    "inverse-proposal-composition"
#define TSEXPR_RSPATH     "proposal-path"
#define TSEXPR_HASTINGS   "hastings-bits"
#define TSEXPR_HTERM      "hastings-term-bits"

#define TSEXPR_TOKEN_PREFIX     "$"
#define TSEXPR_ESC_TOKEN_PREFIX "\\$"   /* escaped version of TSEXPR_TOKEN_PREFIX, for regexps */
#define TSEXPR_BRANCH_INFIX     ":"
#define TSEXPR_WAIT_SUFFIX      "*"

#endif /* TRANSDUCER_KEYWORDS_INCLUDED */
