#ifndef ONTO_KEYWORDS_INCLUDED
#define ONTO_KEYWORDS_INCLUDED

// tags in the Termx data structure
#define TERMX_TERMX    "termx"
#define TERMX_PARAMS         "params"
#define TERMX_MODEL          "model"
#define TERMX_MODEL_SCM      "model-scheme"
#define TERMX_INIT_SCM       "init-scheme"
#define TERMX_TREE_DB_SCM    "tree-db-scheme"
#define TERMX_KNOWLEDGE_SCM  "knowledge-scheme"

// functions provided to Guile Scheme, and their return-value tags
#define GUILE_TERMX_EVIDENCE "termx-evidence"
#define TERMX_LOG_EVIDENCE   "log-evidence"

#define GUILE_TERMX_PREDICT  "termx-predict"
#define TERMX_POSTERIOR      "posterior"

#define GUILE_TERMX_LEARN    "termx-learn"
#define TERMX_TRAINING_LOG   "training-log"
#define TERMX_TRAINING_STEP  "training-step"
#define TERMX_PARAM_COUNTS   "param-counts"
#define TERMX_NEXT_PARAMS    "next-params"

#endif /* ONTO_KEYWORDS_INCLUDED */

