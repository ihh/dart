#ifndef ONTO_KEYWORDS_INCLUDED
#define ONTO_KEYWORDS_INCLUDED

// tags in the Terminatrix data structure
#define TERMINATRIX_TERMINATRIX    "terminatrix"
#define TERMINATRIX_PARAMS         "params"
#define TERMINATRIX_MODEL          "model"
#define TERMINATRIX_MODEL_SCM      "model-scheme"
#define TERMINATRIX_INIT_SCM       "init-scheme"
#define TERMINATRIX_TREE_DB_SCM    "tree-db-scheme"
#define TERMINATRIX_KNOWLEDGE_SCM  "knowledge-scheme"

// functions provided to Guile Scheme, and their return-value tags
#define GUILE_TERMINATRIX_EVIDENCE "terminatrix-evidence"
#define TERMINATRIX_LOG_EVIDENCE   "log-evidence"

#define GUILE_TERMINATRIX_PREDICT  "terminatrix-predict"
#define TERMINATRIX_POSTERIOR      "posterior"

#define GUILE_TERMINATRIX_LEARN    "terminatrix-learn"
#define TERMINATRIX_TRAINING_LOG   "training-log"
#define TERMINATRIX_TRAINING_STEP  "training-step"
#define TERMINATRIX_PARAM_COUNTS   "param-counts"
#define TERMINATRIX_NEXT_PARAMS    "next-params"

#endif /* ONTO_KEYWORDS_INCLUDED */

