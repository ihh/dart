;; Transition matrix
; s_s_s_s (Start)
   s_s_s_s -> s_s_s_il  {(ilp(v))};
   s_s_s_s -> s_s_s_w  {(1-ilp(v))};
; s_s_s_il (Emit: 4)
   s_s_s_il -> s_s_s_il  {(ilp(v))};
   s_s_s_il -> s_s_s_w  {(1-ilp(v)-blp(v))};
   s_s_s_il -> s_s_s_Bi[Ss]  {(blp(v))};
; s_s_s_w (Null)
   s_s_s_w -> s_s_il_w  {(ilp(u))};
   s_s_s_w -> s_s_w_w  {(1-ilp(u))};
; s_s_il_w (Emit: 3)
   s_s_il_w -> s_s_il_w  {(ilp(u))};
   s_s_il_w -> s_s_w_w  {(1-ilp(u)-blp(u))};
   s_s_il_w -> s_s_Bi[Ss]_w  {(blp(u))};
; s_s_w_w (Null)
   s_s_w_w -> s_il_w_w  {(ilp(t))};
   s_s_w_w -> s_w_w_w  {(1-ilp(t))};
; s_il_w_w (Emit: 2)
   s_il_w_w -> s_il_w_w  {(ilp(t))};
   s_il_w_w -> s_w_w_w  {(1-ilp(t)-blp(t))};
   s_il_w_w -> s_Bi[Ss]_w_w  {(blp(t))};
; s_w_w_w (Null)
   s_w_w_w -> il_ml_ml_ml  {(1-blp()-ep()) * (1) * (1) * (1)};
   s_w_w_w -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * (1) * (1) * (1)};
   s_w_w_w -> e_e_e_e  {(ep()) * (1) * (1) * (1)};
; il_ml_ml_s (Null)
   il_ml_ml_s -> il_ml_ml_il  {(ilp(v))};
   il_ml_ml_s -> il_ml_ml_w  {(1-ilp(v))};
; il_ml_ml_il (Emit: 4)
   il_ml_ml_il -> il_ml_ml_il  {(ilp(v))};
   il_ml_ml_il -> il_ml_ml_w  {(1-ilp(v)-blp(v))};
   il_ml_ml_il -> il_ml_ml_Bi[Ss]  {(blp(v))};
; il_ml_ml_ml (Emit: 1,2,3,4)
   il_ml_ml_ml -> il_ml_ml_il  {(ilp(v))};
   il_ml_ml_ml -> il_ml_ml_w  {(1-ilp(v))};
; il_ml_s_w (Null)
   il_ml_s_w -> il_ml_il_w  {(ilp(u))};
   il_ml_s_w -> il_ml_w_w  {(1-ilp(u))};
; il_ml_il_w (Emit: 3)
   il_ml_il_w -> il_ml_il_w  {(ilp(u))};
   il_ml_il_w -> il_ml_w_w  {(1-ilp(u)-blp(u))};
   il_ml_il_w -> il_ml_Bi[Ss]_w  {(blp(u))};
; il_ml_ml_w (Null)
   il_ml_ml_w -> il_ml_il_w  {(ilp(u))};
   il_ml_ml_w -> il_ml_w_w  {(1-ilp(u))};
; il_s_w_w (Null)
   il_s_w_w -> il_il_w_w  {(ilp(t))};
   il_s_w_w -> il_w_w_w  {(1-ilp(t))};
; il_il_w_w (Emit: 2)
   il_il_w_w -> il_il_w_w  {(ilp(t))};
   il_il_w_w -> il_w_w_w  {(1-ilp(t)-blp(t))};
   il_il_w_w -> il_Bi[Ss]_w_w  {(blp(t))};
; il_ml_w_w (Null)
   il_ml_w_w -> il_il_w_w  {(ilp(t))};
   il_ml_w_w -> il_w_w_w  {(1-ilp(t))};
; il_w_w_w (Null)
   il_w_w_w -> il_ml_ml_ml  {(ilp()) * (1) * (1) * (1)};
   il_w_w_w -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * (1) * (1) * (1)};
   il_w_w_w -> e_e_e_e  {(1-ilp()-blp()) * (1) * (1) * (1)};
; e_e_e_S (Start)
   e_e_e_S -> e_e_e_Bi[SS]  {(blp(v))};
   e_e_e_S -> e_e_e_Il  {(1-blp(v)-elp(v))};
   e_e_e_S -> e_e_e_e  {(elp(v))};
; e_e_e_Il (Emit: 4)
   e_e_e_Il -> e_e_e_Bi[SS]  {(blp(v))};
   e_e_e_Il -> e_e_e_Il  {(ilp(v))};
   e_e_e_Il -> e_e_e_e  {(1-ilp(v)-blp(v))};
; e_e_S_e (Start)
   e_e_S_e -> e_e_Bi[SS]_e  {(blp(u))};
   e_e_S_e -> e_e_Il_e  {(1-blp(u)-elp(u))};
   e_e_S_e -> e_e_e_e  {(elp(u))};
; e_e_Il_e (Emit: 3)
   e_e_Il_e -> e_e_Bi[SS]_e  {(blp(u))};
   e_e_Il_e -> e_e_Il_e  {(ilp(u))};
   e_e_Il_e -> e_e_e_e  {(1-ilp(u)-blp(u))};
; e_S_e_e (Start)
   e_S_e_e -> e_Bi[SS]_e_e  {(blp(t))};
   e_S_e_e -> e_Il_e_e  {(1-blp(t)-elp(t))};
   e_S_e_e -> e_e_e_e  {(elp(t))};
; e_Il_e_e (Emit: 2)
   e_Il_e_e -> e_Bi[SS]_e_e  {(blp(t))};
   e_Il_e_e -> e_Il_e_e  {(ilp(t))};
   e_Il_e_e -> e_e_e_e  {(1-ilp(t)-blp(t))};

;; Bifurcations
; s_s_Bi[Ss]_w (Bifurc)
; left:
   s_s_Bi[Ss]_w -> e_e_S_e {(1)};
; center:
   s_s_Bi[Ss]_w -> s_s_s_w {(1)};

; s_s_s_Bi[Ss] (Bifurc)
; left:
   s_s_s_Bi[Ss] -> e_e_e_S {(1)};
; center:
   s_s_s_Bi[Ss] -> s_s_s_s {(1)};

; s_Bi[Ss]_w_w (Bifurc)
; left:
   s_Bi[Ss]_w_w -> e_S_e_e {(1)};
; center:
   s_Bi[Ss]_w_w -> s_s_w_w {(1)};

; il_ml_ml_Bi[Ss] (Bifurc)
; left:
   il_ml_ml_Bi[Ss] -> e_e_e_S {(1)};
; center:
   il_ml_ml_Bi[Ss] -> il_ml_ml_s {(1)};

; il_Bi[Ss]_w_w (Bifurc)
; left:
   il_Bi[Ss]_w_w -> e_S_e_e {(1)};
; center:
   il_Bi[Ss]_w_w -> il_s_w_w {(1)};

; il_ml_Bi[Ss]_w (Bifurc)
; left:
   il_ml_Bi[Ss]_w -> e_e_S_e {(1)};
; center:
   il_ml_Bi[Ss]_w -> il_ml_s_w {(1)};

; Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss] (Bifurc)
; left:
   Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss] -> s_s_s_s {(1) * (1) * (1) * (1)};
; center:
   Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss] -> s_s_s_s {(1) * (1) * (1) * (1)};

; e_e_e_Bi[SS] (Bifurc)
; left:
   e_e_e_Bi[SS] -> e_e_e_S {(1)};
; center:
   e_e_e_Bi[SS] -> e_e_e_S {(1)};

; e_e_Bi[SS]_e (Bifurc)
; left:
   e_e_Bi[SS]_e -> e_e_S_e {(1)};
; center:
   e_e_Bi[SS]_e -> e_e_S_e {(1)};

; e_Bi[SS]_e_e (Bifurc)
; left:
   e_Bi[SS]_e_e -> e_S_e_e {(1)};
; center:
   e_Bi[SS]_e_e -> e_S_e_e {(1)};

; Some statistics on the composed machine(s)

; Singlet machine:
;   4 total states
;   1 deterministic bifurcation states
; Branch machine:
;   10 total states
;   3 deterministic bifurcation states
; Composed machine:
;   34 = 23+10+1 total states
;      4*10*10*10 = 4000 (for comparison)
;   10 deterministic bifurc states
;   23 (non-bifurcation) states in TM
;   60 transitions in TM
;      23^2 = 529 (for comparison)
; Reduced composed machine:
;   30 = 19+10+1 total states
;   10 deterministic bifurc states
;   19 (non-bifurcation) states in reduced TM
;   90 transitions in TM
;      19^2 = 361 (for comparison)
; Fully reduced composed machine:
;   44 = 14+29+1 effective total states
;   29 effective deterministic bifurcation states
;    (10 possibly non-deterministic bifurc states)
;   14 (non-bifurcation) states in TM
;   66 transitions in reduced TM
;      14^2 = 196 (for comparison)
