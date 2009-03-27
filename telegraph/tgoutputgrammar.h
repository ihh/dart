#ifndef TELEGRAPH_GRAMMAR_OUTPUT_ADAPTOR_INCLUDED
#define TELEGRAPH_GRAMMAR_OUTPUT_ADAPTOR_INCLUDED


    /* Output methods */
    void tg_write_grammar (tg_grammar* tg, FILE* file);
    void tg_write_expr (tg_grammar* tg, tg_expr* expr, FILE* file);
    void tg_write_pexpr (tg_grammar* tg, tg_pexpr* pexpr, FILE* file);
    void tg_write_parse (tg_grammar* tg, tg_parse_node* node, FILE* file);
  	 void tg_write_stockholm (tg_grammar* tg, tg_parse_node* node, FILE* file);

#endif /* TELEGRAPH_GRAMMAR_OUTPUT_ADAPTOR_INCLUDED */
