#ifndef CDSPAIRCFG_INCLUDED
#define CDSPAIRCFG_INCLUDED

#include "stemloc/cdspair.h"
#include "scfg/pairpcfg.h"

struct CDS_pair_CFG_state_enum
{
  enum {
    // Dinucleotide-lexicalized states (16 each)
    COD3 = 0,       // type XRYR, lex-emits codon position #3 then goes to COD12
    COD12 = 16,     // type XLRYLR, emits codon positions #1, #2 then goes to End
    INT1COD1 = 32,  // type XLYL, lex-emits codon position #1 then goes to INT1COD23
    INT1COD23 = 48, // type Bifurc, goes to INT and COD23
    COD23 = 64,     // type XLRYLR, emits codon positions #2, #3 then goes to End
    INT2COD3 = 80,  // type XRYR, lex-emits codon position #3 then goes to INT2COD12
    INT2COD12 = 96, // type Bifurc, goes to COD12 and INT

    // Single nucleotide-lexicalized states (4 each)
    // Deletions (X only)
    XCOD3 = 112,      // type XR, lex-emits codon position #3 then goes to XCOD12
    XCOD12 = 116,     // type XLR, emits codon positions #1, #2 then goes to End
    XINT1COD1 = 120,  // type XL, lex-emits codon position #1 then goes to XINT1COD23
    XINT1COD23 = 124, // type Bifurc, goes to XINT and XCOD23
    XCOD23 = 128,     // type XLR, emits codon positions #2, #3 then goes to End
    XINT2COD3 = 132,  // type XR, lex-emits codon position #3 then goes to XINT2COD12
    XINT2COD12 = 136, // type Bifurc, goes to XINT and XCOD12

    // Insertions (Y only)
    YCOD3 = 140,      // type YR, lex-emits codon position #3 then goes to YCOD12
    YCOD12 = 144,     // type YLR, emits codon positions #1, #2 then goes to End
    YINT1COD1 = 148,  // type YL, lex-emits codon position #1 then goes to YINT1COD23
    YINT1COD23 = 152, // type Bifurc, goes to YINT and YCOD23
    YCOD23 = 156,     // type YLR, emits codon positions #2, #3 then goes to End
    YINT2COD3 = 160,  // type YR, lex-emits codon position #3 then goes to YINT2COD12
    YINT2COD12 = 164, // type Bifurc, goes to YINT and YCOD12

    // Unlexicalized states
    // XY codons
    COD = 168,    // type null, codon, goes to COD3 (lexed)
    INTCOD,       // type null, goes to INT0COD or INT1COD1 or INT2COD3
    INT0COD,      // type Bifurc, goes to INT COD
    INT1COD,      // type null, goes to INT1COD1 (lexed)
    INT2COD,      // type null, goes to INT2COD3 (lexed)

    // X codons
    XCOD,         // type null, codon, goes to XCOD3 (lexed)
    XINTCOD,      // type null, goes to XINT0COD or XINT1COD or XINT2COD
    XINT0COD,     // type Bifurc, goes to XINT XCOD
    XINT1COD,     // type null, goes to XINT1COD1 (lexed)
    XINT2COD,     // type null, goes to XINT2COD3 (lexed)

    // Y codons
    YCOD,         // type null, codon, goes to YCOD3 (lexed)
    YINTCOD,      // type null, goes to YINT0COD or YINT1COD1 or YINT2COD3
    YINT0COD,     // type Bifurc, goes to YINT YCOD
    YINT1COD,     // type null, goes to YINT1COD1 (lexed)
    YINT2COD,     // type null, goes to YINT2COD3 (lexed)

    // Introns
    INT,          // type null, goes to XYINT_GG or XINT_GG or YINT_GG
    
    XYINT_GG,     // type XLRYLR, emits first & last intron positions and goes to XYINT_TA
    XYINT_TA,     // type XLRYLR, emits second & penultimate intron positions and goes to XYINTMAT, XYINTINS or XYINTDEL
    XYINTMAT,     // type XLYL, goes to XYINTMAT, XYINTDEL, XYINTINS or End
    XYINTDEL,     // type XL, goes to XYINTMAT, XYINTDEL, XYINTINS or End
    XYINTINS,     // type YL, goes to XYINTMAT, XYINTDEL, XYINTINS or End

    XINT_GG,      // type XLR, emits first & last intron positions and goes to XINT_TA
    XINT_TA,      // type XLR, emits second & penultimate intron positions and goes to XINT
    XINT,         // type XL, goes to XINT or End

    YINT_GG,      // type YLR, emits first & last intron positions and goes to YINT_TA
    YINT_TA,      // type YLR, emits second & penultimate intron positions and goes to YINT
    YINT,         // type YL, goes to YINT or End

    // Exons
    EXMAT,        // type null, goes to CODEXMAT, CODEXINS, CODEXDEL or End
    EXDEL,        // type null, goes to CODEXMAT, CODEXINS, CODEXDEL or End
    EXINS,        // type null, goes to CODEXMAT, CODEXINS, CODEXDEL or End

    CODEXMAT,     // type Bifurc, goes to CODON EXMAT
    CODEXDEL,     // type Bifurc, goes to XCODON EXDEL
    CODEXINS,     // type Bifurc, goes to YCODON EXINS

    CODON,        // type null, goes to COD or INTCOD
    XCODON,       // type null, goes to XCOD or XINTCOD
    YCODON,       // type null, goes to YCOD or YINTCOD
    
    TotalStates   // yields total number of states
  };
};

struct CDS_pair_CFG : CDS_pair_CFG_state_enum, CDS_pair_params, Odds_PCFG
{
  // constructor
  CDS_pair_CFG (PScores& pscore);
};

#endif /* CDSPAIRCFG_INCLUDED */
