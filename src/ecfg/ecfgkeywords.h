#ifndef ECFG_SEXPR_KEYWORDS
#define ECFG_SEXPR_KEYWORDS

#define EG_GRAMMAR        "grammar"
#define EG_META           "meta"
#define EG_TRANSIENT_META "transient-meta"

#define EG_NAME    "name"
#define EG_FROM    "from"
#define EG_TO      "to"
#define EG_PROB    "prob"
#define EG_RATE    "rate"

#define EG_NONTERMINAL "nonterminal"
#define EG_TERMINAL    "terminal"

#define EG_CHAIN            "chain"
#define EG_CHAIN_TERMINAL   EG_TERMINAL
#define EG_CHAIN_POLICY     "update-policy"
#define EG_CHAIN_CLASS      "hidden-class"
#define EG_CHAIN_STATE      "state"
#define EG_CHAIN_INITIAL    "initial"
#define EG_CHAIN_MUTATE     "mutate"

#define EG_POLICY_REV     "rev"
#define EG_POLICY_IRREV   "irrev"
#define EG_POLICY_RIND    "rind"

#define EG_TRANSFORM           "transform"
#define EG_TRANSFORM_ANNOTATE  "annotate"
#define EG_TRANSFORM_ROW       "row"
#define EG_TRANSFORM_COLUMN    "column"
#define EG_TRANSFORM_LABEL     "label"

#define EG_TRANSFORM_MINLEN    "minlen"
#define EG_TRANSFORM_MAXLEN    "maxlen"

#define EG_TRANSFORM_INFIX     "infix"
#define EG_TRANSFORM_PREFIX    "prefix"
#define EG_TRANSFORM_SUFFIX    "suffix"

#define EG_TRANSFORM_SUM_FROM  "sum-from"

#define EG_TRANSFORM_NO_GAPS     "no-gaps"
#define EG_TRANSFORM_STRICT_GAPS "strict-gaps"
#define EG_TRANSFORM_IGNORE_GAPS "gaps-ok"
#define EG_TRANSFORM_GAP_MODEL   "gap-model"
#define EG_TRANSFORM_EXTEND_PROB "extend-prob"
#define EG_TRANSFORM_END_PROB    "end-prob"
#define EG_TRANSFORM_INS_RATE    "insert-rate"
#define EG_TRANSFORM_DEL_RATE    "delete-rate"

#define EG_GFF         "gff"
#define EG_GFF_NONTERM EG_NONTERMINAL
#define EG_GFF_SOURCE  "source"
#define EG_GFF_TYPE    "type"
#define EG_GFF_STRAND  "strand"
#define EG_GFF_FRAME   "frame"

#define EG_UPDATE_RATES "update-rates"
#define EG_UPDATE_RULES "update-rules"

#define EG_PARAMETRIC   "parametric"

// EG_PARAMS and EG_CONST are obsolete, use PK_{PGROUP,RATE,CONST_PGROUP,CONST_RATE} in preference
#define EG_PARAMS       "params"
#define EG_CONST        "const"

#define EG_COUNT        "count"
#define EG_WAIT         "wait"
#define EG_TIME         "time"

#define EG_PSEUDOCOUNTS    "pseudocounts"
#define EG_OBSERVED_COUNTS "observed-counts"
#define EG_CHAIN_COUNTS    "observed-chain-counts"

#define EG_HYBRID_CHAIN      "hybrid-chain"
#define EG_HYBRID_COMPONENTS "components"

// not quite keywords, but important for grammar
#define ECFG_default_start_nonterminal "Start"
#define ECFG_complement_character '~'
#define ECFG_post_emit_character '*'
#define ECFG_deprecated_post_emit_character '\''   /* for backward compatibility */
#define ECFG_default_hidden_class_row_label "hidden_class"

// implicit "#=GS" tags for hybrid chains
#define ECFG_IMPLICIT_GS_TAG_EQUALS   '='    /* prefix for implicit boolean "#=GS" tags defined for each node */
#define ECFG_IMPLICIT_GS_TAG_ANCESTOR ':'    /* prefix for implicit boolean "#=GS" tags defined for each subtree rooted at a particular node */
#define ECFG_IMPLICIT_GS_TAG_NODENAME "?"    /* prefix for implicit "#=GS" tag whose value is the node name */
#define ECFG_IMPLICIT_GS_VALUE_TRUE   "1"    /* "true" value for for implicit boolean "#=GS" tags */
#define ECFG_IMPLICIT_GS_VALUE_FALSE  "0"    /* "false" value for for implicit boolean "#=GS" tags */

// macros
#define EG_FOREACH_TOKEN    "&foreach-token"
#define EG_FOREACH_NODE     "&foreach-node"
#define EG_FOREACH_BRANCH   "&foreach-branch"
#define EG_FOREACH_LEAF     "&foreach-leaf"
#define EG_FOREACH_ANCESTOR "&foreach-ancestor"

#define EG_TOKENS    "#TOKENS"
#define EG_NODES     "#NODES"
#define EG_LEAVES    "#LEAVES"
#define EG_ANCESTORS "#ANCESTORS"
#define EG_BRANCHES  "#BRANCHES"
#define EG_COLUMNS   "#COLUMNS"

#endif /* ECFG_SEXPR_KEYWORDS */
