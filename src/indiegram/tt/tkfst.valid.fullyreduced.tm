;; Fully-reduced transition matrix
; S_S_S_S (Start)
   S_S_S_S -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   S_S_S_S -> S_S_S_IS  {(b2(v))};
   S_S_S_S -> S_S_IS_WS  {(b2(u)) * (1-b2(v))};
   S_S_S_S -> S_IS_WS_WS  {(b2(t)) * (1-b2(u)) * (1-b2(v))};
   S_S_S_S -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   S_S_S_S -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   S_S_S_S -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   S_S_S_S -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   S_S_S_S -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   S_S_S_S -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   S_S_S_S -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   S_S_S_S -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
; L_L_L_L (Start)
   L_L_L_L -> L_L_L_IL  {(b1(v))};
   L_L_L_L -> L_L_IL_WL  {(b1(u)) * (1-b1(v))};
   L_L_L_L -> L_IL_WL_WL  {(b1(t)) * (1-b1(u)) * (1-b1(v))};
   L_L_L_L -> L_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   L_L_L_L -> L_L_Bi[SiL]_WL  {(b1(u)*pS())};
   L_L_L_L -> L_L_L_Bi[SiL]  {(b1(v)*pS())};
   L_L_L_L -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_L -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_L -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_L -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_L -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_L -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_L -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_L -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_L -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_L -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_L -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_L -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_L -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_L -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_L -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_L -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_L -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; S_e_S_S (Start)
   S_e_S_S -> L_e_L_L  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   S_e_S_S -> S_e_S_IS  {(b2(v))};
   S_e_S_S -> S_e_IS_WS  {(b2(u)) * (1-b2(v))};
   S_e_S_S -> IS_e_MS_MS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   S_e_S_S -> IS_e_DS_MS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   S_e_S_S -> IS_e_MS_DS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   S_e_S_S -> IS_e_DS_DS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
; L_e_L_L (Start)
   L_e_L_L -> L_e_L_IL  {(b1(v))};
   L_e_L_L -> L_e_IL_WL  {(b1(u)) * (1-b1(v))};
   L_e_L_L -> L_e_Bi[SiL]_WL  {(b1(u)*pS())};
   L_e_L_L -> L_e_L_Bi[SiL]  {(b1(v)*pS())};
   L_e_L_L -> IL_e_ML_ML  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_e_L_L -> IL_e_DL_ML  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_e_L_L -> IL_e_ML_DL  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_e_L_L -> IL_e_DL_DL  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_e_L_L -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_e_L_L -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_e_L_L -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_e_L_L -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_e_L_L -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; S_S_e_S (Start)
   S_S_e_S -> L_L_e_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   S_S_e_S -> S_S_e_IS  {(b2(v))};
   S_S_e_S -> S_IS_e_WS  {(b2(t)) * (1-b2(v))};
   S_S_e_S -> IS_MS_e_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(v))*(a2(v)))};
   S_S_e_S -> IS_DS_e_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(v))*(a2(v)))};
   S_S_e_S -> IS_MS_e_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(v))*(1-a2(v)))};
   S_S_e_S -> IS_DS_e_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(v))*(1-a2(v)))};
; L_L_e_L (Start)
   L_L_e_L -> L_L_e_IL  {(b1(v))};
   L_L_e_L -> L_IL_e_WL  {(b1(t)) * (1-b1(v))};
   L_L_e_L -> L_Bi[SiL]_e_WL  {(b1(t)*pS())};
   L_L_e_L -> L_L_e_Bi[SiL]  {(b1(v)*pS())};
   L_L_e_L -> IL_ML_e_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   L_L_e_L -> IL_DL_e_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   L_L_e_L -> IL_ML_e_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   L_L_e_L -> IL_DL_e_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   L_L_e_L -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   L_L_e_L -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   L_L_e_L -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   L_L_e_L -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   L_L_e_L -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(v))*(1))};
; S_S_S_IS (Emit: 4)
   S_S_S_IS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   S_S_S_IS -> S_S_S_IS  {(b2(v))};
   S_S_S_IS -> S_S_IS_WS  {(b2(u)) * (1-b2(v))};
   S_S_S_IS -> S_IS_WS_WS  {(b2(t)) * (1-b2(u)) * (1-b2(v))};
   S_S_S_IS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   S_S_S_IS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   S_S_S_IS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   S_S_S_IS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   S_S_S_IS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   S_S_S_IS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   S_S_S_IS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   S_S_S_IS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
; L_L_L_IL (Emit: 4)
   L_L_L_IL -> L_L_L_IL  {(b1(v))};
   L_L_L_IL -> L_L_IL_WL  {(b1(u)) * (1-b1(v))};
   L_L_L_IL -> L_IL_WL_WL  {(b1(t)) * (1-b1(u)) * (1-b1(v))};
   L_L_L_IL -> L_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   L_L_L_IL -> L_L_Bi[SiL]_WL  {(b1(u)*pS())};
   L_L_L_IL -> L_L_L_Bi[SiL]  {(b1(v)*pS())};
   L_L_L_IL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_IL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_IL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_IL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_IL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_IL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_IL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_IL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_L_L_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_L_L_IL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; S_e_S_IS (Emit: 4)
   S_e_S_IS -> L_e_L_L  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   S_e_S_IS -> S_e_S_IS  {(b2(v))};
   S_e_S_IS -> S_e_IS_WS  {(b2(u)) * (1-b2(v))};
   S_e_S_IS -> IS_e_MS_MS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   S_e_S_IS -> IS_e_DS_MS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   S_e_S_IS -> IS_e_MS_DS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   S_e_S_IS -> IS_e_DS_DS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
; L_e_L_IL (Emit: 4)
   L_e_L_IL -> L_e_L_IL  {(b1(v))};
   L_e_L_IL -> L_e_IL_WL  {(b1(u)) * (1-b1(v))};
   L_e_L_IL -> L_e_Bi[SiL]_WL  {(b1(u)*pS())};
   L_e_L_IL -> L_e_L_Bi[SiL]  {(b1(v)*pS())};
   L_e_L_IL -> IL_e_ML_ML  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_e_L_IL -> IL_e_DL_ML  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_e_L_IL -> IL_e_ML_DL  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_e_L_IL -> IL_e_DL_DL  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_e_L_IL -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   L_e_L_IL -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   L_e_L_IL -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_e_L_IL -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   L_e_L_IL -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; L_L_e_IL (Emit: 4)
   L_L_e_IL -> L_L_e_IL  {(b1(v))};
   L_L_e_IL -> L_IL_e_WL  {(b1(t)) * (1-b1(v))};
   L_L_e_IL -> L_Bi[SiL]_e_WL  {(b1(t)*pS())};
   L_L_e_IL -> L_L_e_Bi[SiL]  {(b1(v)*pS())};
   L_L_e_IL -> IL_ML_e_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   L_L_e_IL -> IL_DL_e_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   L_L_e_IL -> IL_ML_e_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   L_L_e_IL -> IL_DL_e_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   L_L_e_IL -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   L_L_e_IL -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   L_L_e_IL -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   L_L_e_IL -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   L_L_e_IL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(v))*(1))};
; L_L_IL_WL (Emit: 3)
   L_L_IL_WL -> L_L_IL_WL  {(b1(u))};
   L_L_IL_WL -> L_IL_WL_WL  {(b1(t)) * (1-b1(u))};
   L_L_IL_WL -> L_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   L_L_IL_WL -> L_L_Bi[SiL]_WL  {(b1(u)*pS())};
   L_L_IL_WL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   L_L_IL_WL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   L_L_IL_WL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   L_L_IL_WL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   L_L_IL_WL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   L_L_IL_WL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   L_L_IL_WL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   L_L_IL_WL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   L_L_IL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   L_L_IL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   L_L_IL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   L_L_IL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   L_L_IL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   L_L_IL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   L_L_IL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   L_L_IL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   L_L_IL_WL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * (1)};
; S_S_IS_WS (Emit: 3)
   S_S_IS_WS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1)*(1))};
   S_S_IS_WS -> S_S_IS_WS  {(b2(u))};
   S_S_IS_WS -> S_IS_WS_WS  {(b2(t)) * (1-b2(u))};
   S_S_IS_WS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * (a2(v))};
   S_S_IS_WS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * (a2(v))};
   S_S_IS_WS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * (a2(v))};
   S_S_IS_WS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * (1-a2(v))};
   S_S_IS_WS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * (a2(v))};
   S_S_IS_WS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * (1-a2(v))};
   S_S_IS_WS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * (1-a2(v))};
   S_S_IS_WS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * (1-a2(v))};
; S_S_S_e (Start)
   S_S_S_e -> L_L_L_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1)))};
   S_S_S_e -> S_S_IS_e  {(b2(u))};
   S_S_S_e -> S_IS_WS_e  {(b2(t)) * (1-b2(u))};
   S_S_S_e -> IS_MS_MS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u)))};
   S_S_S_e -> IS_DS_MS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u)))};
   S_S_S_e -> IS_MS_DS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u)))};
   S_S_S_e -> IS_DS_DS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u)))};
; L_L_L_e (Start)
   L_L_L_e -> L_L_IL_e  {(b1(u))};
   L_L_L_e -> L_IL_WL_e  {(b1(t)) * (1-b1(u))};
   L_L_L_e -> L_Bi[SiL]_WL_e  {(b1(t)*pS())};
   L_L_L_e -> L_L_Bi[SiL]_e  {(b1(u)*pS())};
   L_L_L_e -> IL_ML_ML_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   L_L_L_e -> IL_DL_ML_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   L_L_L_e -> IL_ML_DL_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   L_L_L_e -> IL_DL_DL_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   L_L_L_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   L_L_L_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   L_L_L_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   L_L_L_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   L_L_L_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1))};
; L_IL_WL_WL (Emit: 2)
   L_IL_WL_WL -> L_IL_WL_WL  {(b1(t))};
   L_IL_WL_WL -> L_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   L_IL_WL_WL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v))};
   L_IL_WL_WL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (a1(v))};
   L_IL_WL_WL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (a1(v))};
   L_IL_WL_WL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (1-a1(v))};
   L_IL_WL_WL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (a1(v))};
   L_IL_WL_WL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (1-a1(v))};
   L_IL_WL_WL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (1-a1(v))};
   L_IL_WL_WL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (1-a1(v))};
   L_IL_WL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (1-a1(v))};
   L_IL_WL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (1-a1(v))};
   L_IL_WL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v))};
   L_IL_WL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (a1(v))};
   L_IL_WL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (a1(v))};
   L_IL_WL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (1-a1(v))};
   L_IL_WL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (1-a1(v))};
   L_IL_WL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (a1(v))};
   L_IL_WL_WL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * (1) * (1)};
; S_e_e_S (Start)
   S_e_e_S -> L_e_e_L  {((1-K2())*(1)) * ((1-b2(v))*((1)*(1)))};
   S_e_e_S -> S_e_e_IS  {(b2(v))};
   S_e_e_S -> IS_e_e_MS  {(K2()) * ((1-b2(v))*(a2(v)))};
   S_e_e_S -> IS_e_e_DS  {(K2()) * ((1-b2(v))*(1-a2(v)))};
; L_e_e_L (Start)
   L_e_e_L -> L_e_e_IL  {(b1(v))};
   L_e_e_L -> L_e_e_Bi[SiL]  {(b1(v)*pS())};
   L_e_e_L -> IL_e_e_ML  {(K1()) * ((1-b1(v))*(a1(v)))};
   L_e_e_L -> IL_e_e_DL  {(K1()) * ((1-b1(v))*(1-a1(v)))};
   L_e_e_L -> Bi[SL]_e_e_Bm[SL]  {(K1()*pS()) * ((1-b1(v))*(a1(v)))};
   L_e_e_L -> Bi[SL]_e_e_Bp[eL]  {(K1()*pS()) * ((1-b1(v))*(1-a1(v)))};
   L_e_e_L -> e_e_e_e  {(1-K1()) * ((1-b1(v))*(1))};
; S_S_e_IS (Emit: 4)
   S_S_e_IS -> L_L_e_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   S_S_e_IS -> S_S_e_IS  {(b2(v))};
   S_S_e_IS -> S_IS_e_WS  {(b2(t)) * (1-b2(v))};
   S_S_e_IS -> IS_MS_e_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(v))*(a2(v)))};
   S_S_e_IS -> IS_DS_e_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(v))*(a2(v)))};
   S_S_e_IS -> IS_MS_e_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(v))*(1-a2(v)))};
   S_S_e_IS -> IS_DS_e_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(v))*(1-a2(v)))};
; L_e_IL_WL (Emit: 3)
   L_e_IL_WL -> L_e_IL_WL  {(b1(u))};
   L_e_IL_WL -> L_e_Bi[SiL]_WL  {(b1(u)*pS())};
   L_e_IL_WL -> IL_e_ML_ML  {(K1()) * ((1-b1(u))*(a1(u))) * (a1(v))};
   L_e_IL_WL -> IL_e_DL_ML  {(K1()) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   L_e_IL_WL -> IL_e_ML_DL  {(K1()) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   L_e_IL_WL -> IL_e_DL_DL  {(K1()) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   L_e_IL_WL -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   L_e_IL_WL -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * (a1(v))};
   L_e_IL_WL -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   L_e_IL_WL -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   L_e_IL_WL -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1)) * (1)};
; S_e_IS_WS (Emit: 3)
   S_e_IS_WS -> L_e_L_L  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1))) * ((1)*(1))};
   S_e_IS_WS -> S_e_IS_WS  {(b2(u))};
   S_e_IS_WS -> IS_e_MS_MS  {(K2()) * ((1-b2(u))*(a2(u))) * (a2(v))};
   S_e_IS_WS -> IS_e_DS_MS  {(K2()) * ((1-b2(u))*(1-a2(u))) * (a2(v))};
   S_e_IS_WS -> IS_e_MS_DS  {(K2()) * ((1-b2(u))*(a2(u))) * (1-a2(v))};
   S_e_IS_WS -> IS_e_DS_DS  {(K2()) * ((1-b2(u))*(1-a2(u))) * (1-a2(v))};
; L_IL_e_WL (Emit: 2)
   L_IL_e_WL -> L_IL_e_WL  {(b1(t))};
   L_IL_e_WL -> L_Bi[SiL]_e_WL  {(b1(t)*pS())};
   L_IL_e_WL -> IL_ML_e_ML  {(K1()) * ((1-b1(t))*(a1(t))) * (a1(v))};
   L_IL_e_WL -> IL_DL_e_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * (a1(v))};
   L_IL_e_WL -> IL_ML_e_DL  {(K1()) * ((1-b1(t))*(a1(t))) * (1-a1(v))};
   L_IL_e_WL -> IL_DL_e_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(v))};
   L_IL_e_WL -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(v))};
   L_IL_e_WL -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(v))};
   L_IL_e_WL -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(v))};
   L_IL_e_WL -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(v))};
   L_IL_e_WL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * (1)};
; S_IS_e_WS (Emit: 2)
   S_IS_e_WS -> L_L_e_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1)*(1))};
   S_IS_e_WS -> S_IS_e_WS  {(b2(t))};
   S_IS_e_WS -> IS_MS_e_MS  {(K2()) * ((1-b2(t))*(a2(t))) * (a2(v))};
   S_IS_e_WS -> IS_DS_e_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (a2(v))};
   S_IS_e_WS -> IS_MS_e_DS  {(K2()) * ((1-b2(t))*(a2(t))) * (1-a2(v))};
   S_IS_e_WS -> IS_DS_e_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (1-a2(v))};
; L_e_L_e (Start)
   L_e_L_e -> L_e_IL_e  {(b1(u))};
   L_e_L_e -> L_e_Bi[SiL]_e  {(b1(u)*pS())};
   L_e_L_e -> IL_e_ML_e  {(K1()) * ((1-b1(u))*(a1(u)))};
   L_e_L_e -> IL_e_DL_e  {(K1()) * ((1-b1(u))*(1-a1(u)))};
   L_e_L_e -> Bi[SL]_e_Bp[eL]_e  {(K1()*pS()) * ((1-b1(u))*(1-a1(u)))};
   L_e_L_e -> Bi[SL]_e_Bm[SL]_e  {(K1()*pS()) * ((1-b1(u))*(a1(u)))};
   L_e_L_e -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1))};
; S_e_S_e (Start)
   S_e_S_e -> L_e_L_e  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1)))};
   S_e_S_e -> S_e_IS_e  {(b2(u))};
   S_e_S_e -> IS_e_MS_e  {(K2()) * ((1-b2(u))*(a2(u)))};
   S_e_S_e -> IS_e_DS_e  {(K2()) * ((1-b2(u))*(1-a2(u)))};
; S_S_IS_e (Emit: 3)
   S_S_IS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1)))};
   S_S_IS_e -> S_S_IS_e  {(b2(u))};
   S_S_IS_e -> S_IS_WS_e  {(b2(t)) * (1-b2(u))};
   S_S_IS_e -> IS_MS_MS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u)))};
   S_S_IS_e -> IS_DS_MS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u)))};
   S_S_IS_e -> IS_MS_DS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u)))};
   S_S_IS_e -> IS_DS_DS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u)))};
; L_L_IL_e (Emit: 3)
   L_L_IL_e -> L_L_IL_e  {(b1(u))};
   L_L_IL_e -> L_IL_WL_e  {(b1(t)) * (1-b1(u))};
   L_L_IL_e -> L_Bi[SiL]_WL_e  {(b1(t)*pS())};
   L_L_IL_e -> L_L_Bi[SiL]_e  {(b1(u)*pS())};
   L_L_IL_e -> IL_ML_ML_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   L_L_IL_e -> IL_DL_ML_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   L_L_IL_e -> IL_ML_DL_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   L_L_IL_e -> IL_DL_DL_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   L_L_IL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   L_L_IL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   L_L_IL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   L_L_IL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   L_L_IL_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1))};
; S_IS_WS_WS (Emit: 2)
   S_IS_WS_WS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1)*(1)) * ((1)*(1))};
   S_IS_WS_WS -> S_IS_WS_WS  {(b2(t))};
   S_IS_WS_WS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * (a2(u)) * (a2(v))};
   S_IS_WS_WS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (a2(u)) * (a2(v))};
   S_IS_WS_WS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * (1-a2(u)) * (a2(v))};
   S_IS_WS_WS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * (a2(u)) * (1-a2(v))};
   S_IS_WS_WS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (1-a2(u)) * (a2(v))};
   S_IS_WS_WS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (a2(u)) * (1-a2(v))};
   S_IS_WS_WS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * (1-a2(u)) * (1-a2(v))};
   S_IS_WS_WS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (1-a2(u)) * (1-a2(v))};
; S_e_IS_e (Emit: 3)
   S_e_IS_e -> L_e_L_e  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1)))};
   S_e_IS_e -> S_e_IS_e  {(b2(u))};
   S_e_IS_e -> IS_e_MS_e  {(K2()) * ((1-b2(u))*(a2(u)))};
   S_e_IS_e -> IS_e_DS_e  {(K2()) * ((1-b2(u))*(1-a2(u)))};
; L_e_IL_e (Emit: 3)
   L_e_IL_e -> L_e_IL_e  {(b1(u))};
   L_e_IL_e -> L_e_Bi[SiL]_e  {(b1(u)*pS())};
   L_e_IL_e -> IL_e_ML_e  {(K1()) * ((1-b1(u))*(a1(u)))};
   L_e_IL_e -> IL_e_DL_e  {(K1()) * ((1-b1(u))*(1-a1(u)))};
   L_e_IL_e -> Bi[SL]_e_Bp[eL]_e  {(K1()*pS()) * ((1-b1(u))*(1-a1(u)))};
   L_e_IL_e -> Bi[SL]_e_Bm[SL]_e  {(K1()*pS()) * ((1-b1(u))*(a1(u)))};
   L_e_IL_e -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1))};
; L_IL_WL_e (Emit: 2)
   L_IL_WL_e -> L_IL_WL_e  {(b1(t))};
   L_IL_WL_e -> L_Bi[SiL]_WL_e  {(b1(t)*pS())};
   L_IL_WL_e -> IL_ML_ML_e  {(K1()) * ((1-b1(t))*(a1(t))) * (a1(u))};
   L_IL_WL_e -> IL_DL_ML_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u))};
   L_IL_WL_e -> IL_ML_DL_e  {(K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u))};
   L_IL_WL_e -> IL_DL_DL_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u))};
   L_IL_WL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u))};
   L_IL_WL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u))};
   L_IL_WL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u))};
   L_IL_WL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u))};
   L_IL_WL_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * (1)};
; S_IS_WS_e (Emit: 2)
   S_IS_WS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1)*(1))};
   S_IS_WS_e -> S_IS_WS_e  {(b2(t))};
   S_IS_WS_e -> IS_MS_MS_e  {(K2()) * ((1-b2(t))*(a2(t))) * (a2(u))};
   S_IS_WS_e -> IS_DS_MS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * (a2(u))};
   S_IS_WS_e -> IS_MS_DS_e  {(K2()) * ((1-b2(t))*(a2(t))) * (1-a2(u))};
   S_IS_WS_e -> IS_DS_DS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * (1-a2(u))};
; S_e_e_IS (Emit: 4)
   S_e_e_IS -> L_e_e_L  {((1-K2())*(1)) * ((1-b2(v))*((1)*(1)))};
   S_e_e_IS -> S_e_e_IS  {(b2(v))};
   S_e_e_IS -> IS_e_e_MS  {(K2()) * ((1-b2(v))*(a2(v)))};
   S_e_e_IS -> IS_e_e_DS  {(K2()) * ((1-b2(v))*(1-a2(v)))};
; L_e_e_IL (Emit: 4)
   L_e_e_IL -> L_e_e_IL  {(b1(v))};
   L_e_e_IL -> L_e_e_Bi[SiL]  {(b1(v)*pS())};
   L_e_e_IL -> IL_e_e_ML  {(K1()) * ((1-b1(v))*(a1(v)))};
   L_e_e_IL -> IL_e_e_DL  {(K1()) * ((1-b1(v))*(1-a1(v)))};
   L_e_e_IL -> Bi[SL]_e_e_Bm[SL]  {(K1()*pS()) * ((1-b1(v))*(a1(v)))};
   L_e_e_IL -> Bi[SL]_e_e_Bp[eL]  {(K1()*pS()) * ((1-b1(v))*(1-a1(v)))};
   L_e_e_IL -> e_e_e_e  {(1-K1()) * ((1-b1(v))*(1))};
; L_L_e_e (Start)
   L_L_e_e -> L_IL_e_e  {(b1(t))};
   L_L_e_e -> L_Bi[SiL]_e_e  {(b1(t)*pS())};
   L_L_e_e -> IL_ML_e_e  {(K1()) * ((1-b1(t))*(a1(t)))};
   L_L_e_e -> IL_DL_e_e  {(K1()) * ((1-b1(t))*(1-a1(t)))};
   L_L_e_e -> Bi[SL]_Bp[eL]_e_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t)))};
   L_L_e_e -> Bi[SL]_Bm[SL]_e_e  {(K1()*pS()) * ((1-b1(t))*(a1(t)))};
   L_L_e_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1))};
; S_S_e_e (Start)
   S_S_e_e -> L_L_e_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1)))};
   S_S_e_e -> S_IS_e_e  {(b2(t))};
   S_S_e_e -> IS_MS_e_e  {(K2()) * ((1-b2(t))*(a2(t)))};
   S_S_e_e -> IS_DS_e_e  {(K2()) * ((1-b2(t))*(1-a2(t)))};
; S_IS_e_e (Emit: 2)
   S_IS_e_e -> L_L_e_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1)))};
   S_IS_e_e -> S_IS_e_e  {(b2(t))};
   S_IS_e_e -> IS_MS_e_e  {(K2()) * ((1-b2(t))*(a2(t)))};
   S_IS_e_e -> IS_DS_e_e  {(K2()) * ((1-b2(t))*(1-a2(t)))};
; L_IL_e_e (Emit: 2)
   L_IL_e_e -> L_IL_e_e  {(b1(t))};
   L_IL_e_e -> L_Bi[SiL]_e_e  {(b1(t)*pS())};
   L_IL_e_e -> IL_ML_e_e  {(K1()) * ((1-b1(t))*(a1(t)))};
   L_IL_e_e -> IL_DL_e_e  {(K1()) * ((1-b1(t))*(1-a1(t)))};
   L_IL_e_e -> Bi[SL]_Bp[eL]_e_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t)))};
   L_IL_e_e -> Bi[SL]_Bm[SL]_e_e  {(K1()*pS()) * ((1-b1(t))*(a1(t)))};
   L_IL_e_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1))};
; L_e_e_e (Start)
   L_e_e_e -> IL_e_e_e  {(K1())};
   L_e_e_e -> Bi[SL]_e_e_e  {(K1()*pS())};
   L_e_e_e -> e_e_e_e  {(1-K1())};
; S_e_e_e (Start)
   S_e_e_e -> L_e_e_e  {((1-K2())*(1))};
   S_e_e_e -> IS_e_e_e  {(K2())};
; IS_MS_MS_IS (Emit: 4)
   IS_MS_MS_IS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_MS_MS_IS -> IS_MS_MS_IS  {(b2(v))};
   IS_MS_MS_IS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_MS_IS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_MS_IS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_MS_IS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_MS_IS -> IS_MS_IS_WS  {(b2(u)) * (1-b2(v))};
   IS_MS_MS_IS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_MS_IS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_MS_IS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_MS_IS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_MS_IS -> IS_IS_WS_WS  {(b2(t)) * (1-b2(u)) * (1-b2(v))};
; IL_ML_ML_IL (Emit: 4)
   IL_ML_ML_IL -> IL_ML_ML_IL  {(b1(v))};
   IL_ML_ML_IL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_IL -> IL_ML_IL_WL  {(b1(u)) * (1-b1(v))};
   IL_ML_ML_IL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_IL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_IL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_IL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_IL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_IL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_IL -> IL_IL_WL_WL  {(b1(t)) * (1-b1(u)) * (1-b1(v))};
   IL_ML_ML_IL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_ML_ML_IL -> IL_ML_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_ML_ML_IL -> IL_ML_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_ML_IL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_IL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; IS_MS_DS_IS (Emit: 4)
   IS_MS_DS_IS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-g2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_MS_DS_IS -> IS_MS_DS_IS  {(b2(v))};
   IS_MS_DS_IS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_DS_IS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_DS_IS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_DS_IS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_DS_IS -> IS_MS_IS_WS  {(g2(u)) * (1-b2(v))};
   IS_MS_DS_IS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_DS_IS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_DS_IS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_DS_IS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_DS_IS -> IS_IS_WS_WS  {(b2(t)) * (1-g2(u)) * (1-b2(v))};
; IL_ML_DL_IL (Emit: 4)
   IL_ML_DL_IL -> IL_ML_DL_IL  {(b1(v))};
   IL_ML_DL_IL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_IL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_IL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_IL -> IL_ML_IL_WL  {(g1(u)) * (1-b1(v))};
   IL_ML_DL_IL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_IL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_IL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_IL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_IL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_IL -> IL_IL_WL_WL  {(b1(t)) * (1-g1(u)) * (1-b1(v))};
   IL_ML_DL_IL -> IL_ML_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_ML_DL_IL -> IL_ML_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_DL_IL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_ML_DL_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_IL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-g1(u))*(1)) * ((1-b1(v))*(1))};
; IL_ML_ML_ML (Emit: 1,2,3,4)
   IL_ML_ML_ML -> IL_ML_ML_IL  {(b1(v))};
   IL_ML_ML_ML -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_ML -> IL_ML_IL_WL  {(b1(u)) * (1-b1(v))};
   IL_ML_ML_ML -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_ML -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_ML -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_ML -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_ML -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_ML -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_ML -> IL_IL_WL_WL  {(b1(t)) * (1-b1(u)) * (1-b1(v))};
   IL_ML_ML_ML -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_ML_ML_ML -> IL_ML_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_ML_ML_ML -> IL_ML_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_ML_ML -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_ML -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_ML -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_ML -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_ML -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_ML -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_ML -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_ML -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_ML_ML -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_ML_ML -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; IS_MS_MS_DS (Emit: 1,2,3,4)
   IS_MS_MS_DS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-g2(v))*((1)*(1)))};
   IS_MS_MS_DS -> IS_MS_MS_IS  {(g2(v))};
   IS_MS_MS_DS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_MS_MS_DS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_MS_MS_DS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_MS_MS_DS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_MS_DS -> IS_MS_IS_WS  {(b2(u)) * (1-g2(v))};
   IS_MS_MS_DS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_MS_MS_DS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_MS_DS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_MS_DS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_MS_DS -> IS_IS_WS_WS  {(b2(t)) * (1-b2(u)) * (1-g2(v))};
; IS_MS_IS_WS (Emit: 3)
   IS_MS_IS_WS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1)*(1))};
   IS_MS_IS_WS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * (a2(v))};
   IS_MS_IS_WS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * (a2(v))};
   IS_MS_IS_WS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * (a2(v))};
   IS_MS_IS_WS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * (1-a2(v))};
   IS_MS_IS_WS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * (a2(v))};
   IS_MS_IS_WS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * (1-a2(v))};
   IS_MS_IS_WS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * (1-a2(v))};
   IS_MS_IS_WS -> IS_MS_IS_WS  {(b2(u))};
   IS_MS_IS_WS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * (1-a2(v))};
   IS_MS_IS_WS -> IS_IS_WS_WS  {(b2(t)) * (1-b2(u))};
; IL_ML_IL_e (Emit: 3)
   IL_ML_IL_e -> IL_ML_IL_e  {(b1(u))};
   IL_ML_IL_e -> IL_ML_ML_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_ML_IL_e -> IL_DL_ML_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_ML_IL_e -> IL_ML_DL_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_ML_IL_e -> IL_DL_DL_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_ML_IL_e -> IL_IL_WL_e  {(b1(t)) * (1-b1(u))};
   IL_ML_IL_e -> IL_Bi[SiL]_WL_e  {(b1(t)*pS())};
   IL_ML_IL_e -> IL_ML_Bi[SiL]_e  {(b1(u)*pS())};
   IL_ML_IL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_ML_IL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_ML_IL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_ML_IL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_ML_IL_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1))};
; IS_MS_IS_e (Emit: 3)
   IS_MS_IS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1)))};
   IS_MS_IS_e -> IS_MS_IS_e  {(b2(u))};
   IS_MS_IS_e -> IS_MS_MS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u)))};
   IS_MS_IS_e -> IS_DS_MS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u)))};
   IS_MS_IS_e -> IS_MS_DS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u)))};
   IS_MS_IS_e -> IS_DS_DS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u)))};
   IS_MS_IS_e -> IS_IS_WS_e  {(b2(t)) * (1-b2(u))};
; IS_DS_MS_IS (Emit: 4)
   IS_DS_MS_IS -> L_L_L_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_DS_MS_IS -> IS_DS_MS_IS  {(b2(v))};
   IS_DS_MS_IS -> IS_MS_MS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_MS_IS -> IS_DS_MS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_MS_IS -> IS_MS_DS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_MS_IS -> IS_MS_MS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_MS_IS -> IS_DS_DS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_MS_IS -> IS_DS_MS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_MS_IS -> IS_MS_DS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_MS_IS -> IS_DS_DS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_MS_IS -> IS_DS_IS_WS  {(b2(u)) * (1-b2(v))};
   IS_DS_MS_IS -> IS_IS_WS_WS  {(g2(t)) * (1-b2(u)) * (1-b2(v))};
; IL_DL_ML_IL (Emit: 4)
   IL_DL_ML_IL -> IL_DL_ML_IL  {(b1(v))};
   IL_DL_ML_IL -> IL_ML_ML_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_IL -> IL_DL_ML_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_IL -> IL_ML_DL_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_IL -> IL_ML_ML_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_IL -> IL_DL_DL_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_IL -> IL_DL_ML_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_IL -> IL_ML_DL_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_IL -> IL_DL_DL_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_IL -> IL_DL_IL_WL  {(b1(u)) * (1-b1(v))};
   IL_DL_ML_IL -> IL_IL_WL_WL  {(g1(t)) * (1-b1(u)) * (1-b1(v))};
   IL_DL_ML_IL -> IL_DL_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_ML_IL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_DL_ML_IL -> IL_DL_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_DL_ML_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_IL -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; IS_e_MS_IS (Emit: 4)
   IS_e_MS_IS -> L_e_L_L  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_e_MS_IS -> IS_e_MS_IS  {(b2(v))};
   IS_e_MS_IS -> IS_e_MS_MS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_e_MS_IS -> IS_e_DS_MS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_e_MS_IS -> IS_e_MS_DS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_e_MS_IS -> IS_e_DS_DS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_e_MS_IS -> IS_e_IS_WS  {(b2(u)) * (1-b2(v))};
; IL_e_ML_IL (Emit: 4)
   IL_e_ML_IL -> IL_e_ML_IL  {(b1(v))};
   IL_e_ML_IL -> IL_e_ML_ML  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_ML_IL -> IL_e_DL_ML  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_ML_IL -> IL_e_ML_DL  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_ML_IL -> IL_e_DL_DL  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_ML_IL -> IL_e_IL_WL  {(b1(u)) * (1-b1(v))};
   IL_e_ML_IL -> IL_e_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_e_ML_IL -> IL_e_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_e_ML_IL -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_ML_IL -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_ML_IL -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_ML_IL -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_ML_IL -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; IS_MS_MS_MS (Emit: 1,2,3,4)
   IS_MS_MS_MS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_MS_MS_MS -> IS_MS_MS_IS  {(b2(v))};
   IS_MS_MS_MS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_MS_MS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_MS_MS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_MS_MS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_MS_MS -> IS_MS_IS_WS  {(b2(u)) * (1-b2(v))};
   IS_MS_MS_MS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_MS_MS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_MS_MS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_MS_MS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_MS_MS -> IS_IS_WS_WS  {(b2(t)) * (1-b2(u)) * (1-b2(v))};
; IL_DL_DL_IL (Emit: 4)
   IL_DL_DL_IL -> IL_ML_ML_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_IL -> IL_DL_ML_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_IL -> IL_ML_DL_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_IL -> IL_ML_ML_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_IL -> IL_DL_DL_IL  {(b1(v))};
   IL_DL_DL_IL -> IL_DL_DL_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_IL -> IL_DL_ML_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_IL -> IL_ML_DL_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_IL -> IL_DL_DL_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_IL -> IL_DL_IL_WL  {(g1(u)) * (1-b1(v))};
   IL_DL_DL_IL -> IL_IL_WL_WL  {(g1(t)) * (1-g1(u)) * (1-b1(v))};
   IL_DL_DL_IL -> IL_DL_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_DL_DL_IL -> IL_DL_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_DL_IL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_DL_DL_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_IL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_IL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_IL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_IL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_IL -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-g1(u))*(1)) * ((1-b1(v))*(1))};
; IS_DS_DS_IS (Emit: 4)
   IS_DS_DS_IS -> L_L_L_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-g2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_DS_DS_IS -> IS_MS_MS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_DS_IS -> IS_DS_DS_IS  {(b2(v))};
   IS_DS_DS_IS -> IS_DS_MS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_DS_IS -> IS_MS_DS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_DS_IS -> IS_MS_MS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_DS_IS -> IS_DS_DS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_DS_IS -> IS_DS_MS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_DS_IS -> IS_MS_DS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_DS_IS -> IS_DS_DS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_DS_IS -> IS_DS_IS_WS  {(g2(u)) * (1-b2(v))};
   IS_DS_DS_IS -> IS_IS_WS_WS  {(g2(t)) * (1-g2(u)) * (1-b2(v))};
; IL_ML_e_IL (Emit: 4)
   IL_ML_e_IL -> IL_ML_e_IL  {(b1(v))};
   IL_ML_e_IL -> IL_ML_e_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_ML_e_IL -> IL_DL_e_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_ML_e_IL -> IL_ML_e_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_e_IL -> IL_DL_e_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_e_IL -> IL_IL_e_WL  {(b1(t)) * (1-b1(v))};
   IL_ML_e_IL -> IL_Bi[SiL]_e_WL  {(b1(t)*pS())};
   IL_ML_e_IL -> IL_ML_e_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_e_IL -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_ML_e_IL -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_e_IL -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_e_IL -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_ML_e_IL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(v))*(1))};
; IS_MS_e_IS (Emit: 4)
   IS_MS_e_IS -> L_L_e_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_MS_e_IS -> IS_MS_e_IS  {(b2(v))};
   IS_MS_e_IS -> IS_MS_e_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(v))*(a2(v)))};
   IS_MS_e_IS -> IS_DS_e_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(v))*(a2(v)))};
   IS_MS_e_IS -> IS_MS_e_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_e_IS -> IS_DS_e_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_e_IS -> IS_IS_e_WS  {(b2(t)) * (1-b2(v))};
; IL_DL_e_IL (Emit: 4)
   IL_DL_e_IL -> IL_DL_e_IL  {(b1(v))};
   IL_DL_e_IL -> IL_ML_e_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_DL_e_IL -> IL_DL_e_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_DL_e_IL -> IL_ML_e_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_e_IL -> IL_DL_e_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_e_IL -> IL_IL_e_WL  {(g1(t)) * (1-b1(v))};
   IL_DL_e_IL -> IL_Bi[SiL]_e_WL  {(b1(t)*pS())};
   IL_DL_e_IL -> IL_DL_e_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_e_IL -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_DL_e_IL -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_e_IL -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_e_IL -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_DL_e_IL -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-b1(v))*(1))};
; IS_DS_MS_MS (Emit: 1,2,3,4)
   IS_DS_MS_MS -> L_L_L_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_DS_MS_MS -> IS_DS_MS_IS  {(b2(v))};
   IS_DS_MS_MS -> IS_MS_MS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_MS_MS -> IS_DS_MS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_MS_MS -> IS_MS_DS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_MS_MS -> IS_MS_MS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_MS_MS -> IS_DS_DS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_MS_MS -> IS_DS_MS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_MS_MS -> IS_MS_DS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_MS_MS -> IS_DS_DS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_MS_MS -> IS_DS_IS_WS  {(b2(u)) * (1-b2(v))};
   IS_DS_MS_MS -> IS_IS_WS_WS  {(g2(t)) * (1-b2(u)) * (1-b2(v))};
; IL_DL_ML_ML (Emit: 1,2,3,4)
   IL_DL_ML_ML -> IL_DL_ML_IL  {(b1(v))};
   IL_DL_ML_ML -> IL_ML_ML_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_ML -> IL_DL_ML_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_ML -> IL_ML_DL_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_ML -> IL_ML_ML_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_ML -> IL_DL_DL_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_ML -> IL_DL_ML_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_ML -> IL_ML_DL_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_ML -> IL_DL_DL_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_ML -> IL_DL_IL_WL  {(b1(u)) * (1-b1(v))};
   IL_DL_ML_ML -> IL_IL_WL_WL  {(g1(t)) * (1-b1(u)) * (1-b1(v))};
   IL_DL_ML_ML -> IL_DL_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_ML_ML -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_DL_ML_ML -> IL_DL_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_DL_ML_ML -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_ML -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_ML -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_ML -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_ML -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_ML -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_ML -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_ML_ML -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_ML_ML -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; IS_MS_DS_MS (Emit: 1,2,3,4)
   IS_MS_DS_MS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-g2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_MS_DS_MS -> IS_MS_DS_IS  {(b2(v))};
   IS_MS_DS_MS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_DS_MS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_DS_MS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_DS_MS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_DS_MS -> IS_MS_IS_WS  {(g2(u)) * (1-b2(v))};
   IS_MS_DS_MS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_MS_DS_MS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_DS_MS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_DS_MS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_DS_MS -> IS_IS_WS_WS  {(b2(t)) * (1-g2(u)) * (1-b2(v))};
; IL_ML_ML_DL (Emit: 1,2,3,4)
   IL_ML_ML_DL -> IL_ML_ML_IL  {(g1(v))};
   IL_ML_ML_DL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_ML_DL -> IL_ML_IL_WL  {(b1(u)) * (1-g1(v))};
   IL_ML_ML_DL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_ML_DL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_ML_DL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_ML_DL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_ML_DL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_ML_DL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_ML_DL -> IL_IL_WL_WL  {(b1(t)) * (1-b1(u)) * (1-g1(v))};
   IL_ML_ML_DL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_ML_ML_DL -> IL_ML_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_ML_ML_DL -> IL_ML_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_ML_DL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_ML_DL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_ML_DL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_ML_DL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_ML_DL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_ML_DL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_ML_DL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_ML_DL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_ML_DL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_ML_DL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * ((1-g1(v))*(1))};
; IS_DS_MS_DS (Emit: 1,2,3,4)
   IS_DS_MS_DS -> L_L_L_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1-g2(v))*((1)*(1)))};
   IS_DS_MS_DS -> IS_DS_MS_IS  {(g2(v))};
   IS_DS_MS_DS -> IS_MS_MS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_DS_MS_DS -> IS_DS_MS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_DS_MS_DS -> IS_MS_DS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_DS_MS_DS -> IS_MS_MS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_MS_DS -> IS_DS_DS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_DS_MS_DS -> IS_DS_MS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_MS_DS -> IS_MS_DS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_MS_DS -> IS_DS_DS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_MS_DS -> IS_DS_IS_WS  {(b2(u)) * (1-g2(v))};
   IS_DS_MS_DS -> IS_IS_WS_WS  {(g2(t)) * (1-b2(u)) * (1-g2(v))};
; IL_e_DL_IL (Emit: 4)
   IL_e_DL_IL -> IL_e_DL_IL  {(b1(v))};
   IL_e_DL_IL -> IL_e_ML_ML  {(K1()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_DL_IL -> IL_e_DL_ML  {(K1()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_DL_IL -> IL_e_ML_DL  {(K1()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_DL_IL -> IL_e_DL_DL  {(K1()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_DL_IL -> IL_e_IL_WL  {(g1(u)) * (1-b1(v))};
   IL_e_DL_IL -> IL_e_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_e_DL_IL -> IL_e_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_e_DL_IL -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_DL_IL -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_DL_IL -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_DL_IL -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_DL_IL -> e_e_e_e  {(1-K1()) * ((1-g1(u))*(1)) * ((1-b1(v))*(1))};
; IS_e_DS_IS (Emit: 4)
   IS_e_DS_IS -> L_e_L_L  {((1-K2())*(1)) * ((1-g2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_e_DS_IS -> IS_e_DS_IS  {(b2(v))};
   IS_e_DS_IS -> IS_e_MS_MS  {(K2()) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_e_DS_IS -> IS_e_DS_MS  {(K2()) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_e_DS_IS -> IS_e_MS_DS  {(K2()) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_e_DS_IS -> IS_e_DS_DS  {(K2()) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_e_DS_IS -> IS_e_IS_WS  {(g2(u)) * (1-b2(v))};
; IL_ML_DL_ML (Emit: 1,2,3,4)
   IL_ML_DL_ML -> IL_ML_DL_IL  {(b1(v))};
   IL_ML_DL_ML -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_ML -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_ML -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_ML -> IL_ML_IL_WL  {(g1(u)) * (1-b1(v))};
   IL_ML_DL_ML -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_ML -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_ML -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_ML -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_ML -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_ML -> IL_IL_WL_WL  {(b1(t)) * (1-g1(u)) * (1-b1(v))};
   IL_ML_DL_ML -> IL_ML_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_ML_DL_ML -> IL_ML_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_DL_ML -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_ML_DL_ML -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_ML -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_ML -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_ML -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_ML -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_ML -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_ML -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_DL_ML -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_ML_DL_ML -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-g1(u))*(1)) * ((1-b1(v))*(1))};
; IS_DS_e_IS (Emit: 4)
   IS_DS_e_IS -> L_L_e_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_DS_e_IS -> IS_DS_e_IS  {(b2(v))};
   IS_DS_e_IS -> IS_MS_e_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(v))*(a2(v)))};
   IS_DS_e_IS -> IS_DS_e_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(v))*(a2(v)))};
   IS_DS_e_IS -> IS_MS_e_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_e_IS -> IS_DS_e_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_e_IS -> IS_IS_e_WS  {(g2(t)) * (1-b2(v))};
; IL_e_ML_ML (Emit: 1,3,4)
   IL_e_ML_ML -> IL_e_ML_IL  {(b1(v))};
   IL_e_ML_ML -> IL_e_ML_ML  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_ML_ML -> IL_e_DL_ML  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_ML_ML -> IL_e_ML_DL  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_ML_ML -> IL_e_DL_DL  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_ML_ML -> IL_e_IL_WL  {(b1(u)) * (1-b1(v))};
   IL_e_ML_ML -> IL_e_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_e_ML_ML -> IL_e_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_e_ML_ML -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_ML_ML -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_ML_ML -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_ML_ML -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_ML_ML -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1)) * ((1-b1(v))*(1))};
; IL_DL_DL_ML (Emit: 1,2,3,4)
   IL_DL_DL_ML -> IL_ML_ML_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_ML -> IL_DL_ML_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_ML -> IL_ML_DL_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_ML -> IL_ML_ML_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_ML -> IL_DL_DL_IL  {(b1(v))};
   IL_DL_DL_ML -> IL_DL_DL_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_ML -> IL_DL_ML_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_ML -> IL_ML_DL_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_ML -> IL_DL_DL_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_ML -> IL_DL_IL_WL  {(g1(u)) * (1-b1(v))};
   IL_DL_DL_ML -> IL_IL_WL_WL  {(g1(t)) * (1-g1(u)) * (1-b1(v))};
   IL_DL_DL_ML -> IL_DL_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_DL_DL_ML -> IL_DL_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_DL_ML -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_DL_DL_ML -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_ML -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_ML -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_ML -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_ML -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_ML -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_ML -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_DL_ML -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_DL_DL_ML -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-g1(u))*(1)) * ((1-b1(v))*(1))};
; IL_ML_e_ML (Emit: 1,2,4)
   IL_ML_e_ML -> IL_ML_e_IL  {(b1(v))};
   IL_ML_e_ML -> IL_ML_e_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_ML_e_ML -> IL_DL_e_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_ML_e_ML -> IL_ML_e_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_e_ML -> IL_DL_e_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_e_ML -> IL_IL_e_WL  {(b1(t)) * (1-b1(v))};
   IL_ML_e_ML -> IL_Bi[SiL]_e_WL  {(b1(t)*pS())};
   IL_ML_e_ML -> IL_ML_e_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_e_ML -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_ML_e_ML -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_e_ML -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_ML_e_ML -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_ML_e_ML -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(v))*(1))};
; IL_DL_ML_DL (Emit: 1,2,3,4)
   IL_DL_ML_DL -> IL_DL_ML_IL  {(g1(v))};
   IL_DL_ML_DL -> IL_ML_ML_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_ML_DL -> IL_DL_ML_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_ML_DL -> IL_ML_DL_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_ML_DL -> IL_ML_ML_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_ML_DL -> IL_DL_DL_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_ML_DL -> IL_DL_ML_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_ML_DL -> IL_ML_DL_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_ML_DL -> IL_DL_DL_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_ML_DL -> IL_DL_IL_WL  {(b1(u)) * (1-g1(v))};
   IL_DL_ML_DL -> IL_IL_WL_WL  {(g1(t)) * (1-b1(u)) * (1-g1(v))};
   IL_DL_ML_DL -> IL_DL_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_ML_DL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_DL_ML_DL -> IL_DL_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_DL_ML_DL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_ML_DL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_ML_DL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_ML_DL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_ML_DL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_ML_DL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_ML_DL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_ML_DL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_ML_DL -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1)) * ((1-g1(v))*(1))};
; IS_e_DS_MS (Emit: 1,3,4)
   IS_e_DS_MS -> L_e_L_L  {((1-K2())*(1)) * ((1-g2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_e_DS_MS -> IS_e_DS_IS  {(b2(v))};
   IS_e_DS_MS -> IS_e_MS_MS  {(K2()) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_e_DS_MS -> IS_e_DS_MS  {(K2()) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_e_DS_MS -> IS_e_MS_DS  {(K2()) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_e_DS_MS -> IS_e_DS_DS  {(K2()) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_e_DS_MS -> IS_e_IS_WS  {(g2(u)) * (1-b2(v))};
; IS_e_MS_DS (Emit: 1,3,4)
   IS_e_MS_DS -> L_e_L_L  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1))) * ((1-g2(v))*((1)*(1)))};
   IS_e_MS_DS -> IS_e_MS_IS  {(g2(v))};
   IS_e_MS_DS -> IS_e_MS_MS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_e_MS_DS -> IS_e_DS_MS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_e_MS_DS -> IS_e_MS_DS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_e_MS_DS -> IS_e_DS_DS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_e_MS_DS -> IS_e_IS_WS  {(b2(u)) * (1-g2(v))};
; IL_ML_DL_DL (Emit: 1,2,3,4)
   IL_ML_DL_DL -> IL_ML_DL_IL  {(g1(v))};
   IL_ML_DL_DL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_DL_DL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_DL_DL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_DL_DL -> IL_ML_IL_WL  {(g1(u)) * (1-g1(v))};
   IL_ML_DL_DL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_DL_DL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_DL_DL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_DL_DL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_DL_DL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_DL_DL -> IL_IL_WL_WL  {(b1(t)) * (1-g1(u)) * (1-g1(v))};
   IL_ML_DL_DL -> IL_ML_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_ML_DL_DL -> IL_ML_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_DL_DL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_ML_DL_DL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_DL_DL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_DL_DL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_DL_DL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_DL_DL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_DL_DL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_DL_DL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_DL_DL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_ML_DL_DL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-g1(u))*(1)) * ((1-g1(v))*(1))};
; IS_MS_DS_DS (Emit: 1,2,3,4)
   IS_MS_DS_DS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-g2(u))*((1)*(1))) * ((1-g2(v))*((1)*(1)))};
   IS_MS_DS_DS -> IS_MS_DS_IS  {(g2(v))};
   IS_MS_DS_DS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_MS_DS_DS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_MS_DS_DS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_MS_DS_DS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_DS_DS -> IS_MS_IS_WS  {(g2(u)) * (1-g2(v))};
   IS_MS_DS_DS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_MS_DS_DS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_DS_DS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_DS_DS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_DS_DS -> IS_IS_WS_WS  {(b2(t)) * (1-g2(u)) * (1-g2(v))};
; IL_ML_ML_e (Emit: 1,2,3)
   IL_ML_ML_e -> IL_ML_IL_e  {(b1(u))};
   IL_ML_ML_e -> IL_ML_ML_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_ML_ML_e -> IL_DL_ML_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_ML_ML_e -> IL_ML_DL_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_ML_ML_e -> IL_DL_DL_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_ML_ML_e -> IL_IL_WL_e  {(b1(t)) * (1-b1(u))};
   IL_ML_ML_e -> IL_Bi[SiL]_WL_e  {(b1(t)*pS())};
   IL_ML_ML_e -> IL_ML_Bi[SiL]_e  {(b1(u)*pS())};
   IL_ML_ML_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_ML_ML_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_ML_ML_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_ML_ML_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_ML_ML_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1))};
; IS_MS_MS_e (Emit: 1,2,3)
   IS_MS_MS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1)))};
   IS_MS_MS_e -> IS_MS_IS_e  {(b2(u))};
   IS_MS_MS_e -> IS_MS_MS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(a2(u)))};
   IS_MS_MS_e -> IS_DS_MS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u)))};
   IS_MS_MS_e -> IS_MS_DS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u)))};
   IS_MS_MS_e -> IS_DS_DS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u)))};
   IS_MS_MS_e -> IS_IS_WS_e  {(b2(t)) * (1-b2(u))};
; IS_DS_DS_MS (Emit: 1,2,3,4)
   IS_DS_DS_MS -> L_L_L_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-g2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_DS_DS_MS -> IS_MS_MS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_DS_MS -> IS_DS_DS_IS  {(b2(v))};
   IS_DS_DS_MS -> IS_DS_MS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_DS_MS -> IS_MS_DS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_DS_MS -> IS_MS_MS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_DS_MS -> IS_DS_DS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_DS_DS_MS -> IS_DS_MS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_DS_MS -> IS_MS_DS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_DS_MS -> IS_DS_DS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_DS_MS -> IS_DS_IS_WS  {(g2(u)) * (1-b2(v))};
   IS_DS_DS_MS -> IS_IS_WS_WS  {(g2(t)) * (1-g2(u)) * (1-b2(v))};
; IS_DS_e_MS (Emit: 1,2,4)
   IS_DS_e_MS -> L_L_e_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_DS_e_MS -> IS_DS_e_IS  {(b2(v))};
   IS_DS_e_MS -> IS_MS_e_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(v))*(a2(v)))};
   IS_DS_e_MS -> IS_DS_e_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(v))*(a2(v)))};
   IS_DS_e_MS -> IS_MS_e_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_e_MS -> IS_DS_e_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(v))*(1-a2(v)))};
   IS_DS_e_MS -> IS_IS_e_WS  {(g2(t)) * (1-b2(v))};
; IL_DL_DL_DL (Emit: 1,2,3,4)
   IL_DL_DL_DL -> IL_ML_ML_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_DL_DL -> IL_DL_ML_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_DL_DL -> IL_ML_DL_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_DL_DL -> IL_ML_ML_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_DL_DL -> IL_DL_DL_IL  {(g1(v))};
   IL_DL_DL_DL -> IL_DL_DL_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_DL_DL -> IL_DL_ML_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_DL_DL -> IL_ML_DL_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_DL_DL -> IL_DL_DL_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_DL_DL -> IL_DL_IL_WL  {(g1(u)) * (1-g1(v))};
   IL_DL_DL_DL -> IL_IL_WL_WL  {(g1(t)) * (1-g1(u)) * (1-g1(v))};
   IL_DL_DL_DL -> IL_DL_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_DL_DL_DL -> IL_DL_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_DL_DL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_DL_DL_DL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_DL_DL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_DL_DL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_DL_DL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_DL_DL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_DL_DL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_DL_DL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_DL_DL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_DL_DL_DL -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-g1(u))*(1)) * ((1-g1(v))*(1))};
; IS_DS_DS_DS (Emit: 1,2,3,4)
   IS_DS_DS_DS -> L_L_L_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-g2(u))*((1)*(1))) * ((1-g2(v))*((1)*(1)))};
   IS_DS_DS_DS -> IS_MS_MS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_DS_DS_DS -> IS_DS_DS_IS  {(g2(v))};
   IS_DS_DS_DS -> IS_DS_MS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_DS_DS_DS -> IS_MS_DS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_DS_DS_DS -> IS_MS_MS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_DS_DS -> IS_DS_DS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_DS_DS_DS -> IS_DS_MS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_DS_DS -> IS_MS_DS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_DS_DS -> IS_DS_DS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_DS_DS -> IS_DS_IS_WS  {(g2(u)) * (1-g2(v))};
   IS_DS_DS_DS -> IS_IS_WS_WS  {(g2(t)) * (1-g2(u)) * (1-g2(v))};
; IL_DL_IL_WL (Emit: 3)
   IL_DL_IL_WL -> IL_ML_ML_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_DL_IL_WL -> IL_DL_ML_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_DL_IL_WL -> IL_ML_DL_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_DL_IL_WL -> IL_ML_ML_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_DL_IL_WL -> IL_DL_DL_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_DL_IL_WL -> IL_DL_ML_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_DL_IL_WL -> IL_ML_DL_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_DL_IL_WL -> IL_DL_DL_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_DL_IL_WL -> IL_DL_IL_WL  {(b1(u))};
   IL_DL_IL_WL -> IL_IL_WL_WL  {(g1(t)) * (1-b1(u))};
   IL_DL_IL_WL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_DL_IL_WL -> IL_DL_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_DL_IL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_DL_IL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_DL_IL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_DL_IL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_DL_IL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_DL_IL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_DL_IL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_DL_IL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_DL_IL_WL -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1)) * (1)};
; IL_IL_WL_WL (Emit: 2)
   IL_IL_WL_WL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v))};
   IL_IL_WL_WL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (a1(v))};
   IL_IL_WL_WL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (a1(v))};
   IL_IL_WL_WL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (1-a1(v))};
   IL_IL_WL_WL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (a1(v))};
   IL_IL_WL_WL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (1-a1(v))};
   IL_IL_WL_WL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (1-a1(v))};
   IL_IL_WL_WL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (1-a1(v))};
   IL_IL_WL_WL -> IL_IL_WL_WL  {(b1(t))};
   IL_IL_WL_WL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_IL_WL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (1-a1(v))};
   IL_IL_WL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (1-a1(v))};
   IL_IL_WL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v))};
   IL_IL_WL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (a1(v))};
   IL_IL_WL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (a1(v))};
   IL_IL_WL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (1-a1(v))};
   IL_IL_WL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (1-a1(v))};
   IL_IL_WL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (a1(v))};
   IL_IL_WL_WL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * (1) * (1)};
; IL_DL_ML_e (Emit: 1,2,3)
   IL_DL_ML_e -> IL_DL_IL_e  {(b1(u))};
   IL_DL_ML_e -> IL_ML_ML_e  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_DL_ML_e -> IL_DL_ML_e  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_DL_ML_e -> IL_ML_DL_e  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_DL_ML_e -> IL_DL_DL_e  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_DL_ML_e -> IL_IL_WL_e  {(g1(t)) * (1-b1(u))};
   IL_DL_ML_e -> IL_Bi[SiL]_WL_e  {(b1(t)*pS())};
   IL_DL_ML_e -> IL_DL_Bi[SiL]_e  {(b1(u)*pS())};
   IL_DL_ML_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_DL_ML_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_DL_ML_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_DL_ML_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_DL_ML_e -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1))};
; IL_e_DL_ML (Emit: 1,3,4)
   IL_e_DL_ML -> IL_e_DL_IL  {(b1(v))};
   IL_e_DL_ML -> IL_e_ML_ML  {(K1()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_DL_ML -> IL_e_DL_ML  {(K1()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_DL_ML -> IL_e_ML_DL  {(K1()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_DL_ML -> IL_e_DL_DL  {(K1()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_DL_ML -> IL_e_IL_WL  {(g1(u)) * (1-b1(v))};
   IL_e_DL_ML -> IL_e_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_e_DL_ML -> IL_e_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_e_DL_ML -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_DL_ML -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v)))};
   IL_e_DL_ML -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_DL_ML -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v)))};
   IL_e_DL_ML -> e_e_e_e  {(1-K1()) * ((1-g1(u))*(1)) * ((1-b1(v))*(1))};
; IS_e_e_MS (Emit: 1,4)
   IS_e_e_MS -> L_e_e_L  {((1-K2())*(1)) * ((1-b2(v))*((1)*(1)))};
   IS_e_e_MS -> IS_e_e_IS  {(b2(v))};
   IS_e_e_MS -> IS_e_e_MS  {(K2()) * ((1-b2(v))*(a2(v)))};
   IS_e_e_MS -> IS_e_e_DS  {(K2()) * ((1-b2(v))*(1-a2(v)))};
; IL_e_e_ML (Emit: 1,4)
   IL_e_e_ML -> IL_e_e_IL  {(b1(v))};
   IL_e_e_ML -> IL_e_e_ML  {(K1()) * ((1-b1(v))*(a1(v)))};
   IL_e_e_ML -> IL_e_e_DL  {(K1()) * ((1-b1(v))*(1-a1(v)))};
   IL_e_e_ML -> IL_e_e_Bi[SiL]  {(b1(v)*pS())};
   IL_e_e_ML -> Bi[SL]_e_e_Bm[SL]  {(K1()*pS()) * ((1-b1(v))*(a1(v)))};
   IL_e_e_ML -> Bi[SL]_e_e_Bp[eL]  {(K1()*pS()) * ((1-b1(v))*(1-a1(v)))};
   IL_e_e_ML -> e_e_e_e  {(1-K1()) * ((1-b1(v))*(1))};
; IL_e_DL_DL (Emit: 1,3,4)
   IL_e_DL_DL -> IL_e_DL_IL  {(g1(v))};
   IL_e_DL_DL -> IL_e_ML_ML  {(K1()) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_e_DL_DL -> IL_e_DL_ML  {(K1()) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_e_DL_DL -> IL_e_ML_DL  {(K1()) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_e_DL_DL -> IL_e_DL_DL  {(K1()) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_e_DL_DL -> IL_e_IL_WL  {(g1(u)) * (1-g1(v))};
   IL_e_DL_DL -> IL_e_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_e_DL_DL -> IL_e_DL_Bi[SiL]  {(b1(v)*pS())};
   IL_e_DL_DL -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_e_DL_DL -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_e_DL_DL -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-g1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_e_DL_DL -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-g1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_e_DL_DL -> e_e_e_e  {(1-K1()) * ((1-g1(u))*(1)) * ((1-g1(v))*(1))};
; IL_ML_e_DL (Emit: 1,2,4)
   IL_ML_e_DL -> IL_ML_e_IL  {(g1(v))};
   IL_ML_e_DL -> IL_ML_e_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(v))*(a1(v)))};
   IL_ML_e_DL -> IL_DL_e_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(v))*(a1(v)))};
   IL_ML_e_DL -> IL_ML_e_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_e_DL -> IL_DL_e_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_e_DL -> IL_IL_e_WL  {(b1(t)) * (1-g1(v))};
   IL_ML_e_DL -> IL_Bi[SiL]_e_WL  {(b1(t)*pS())};
   IL_ML_e_DL -> IL_ML_e_Bi[SiL]  {(b1(v)*pS())};
   IL_ML_e_DL -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(v))*(a1(v)))};
   IL_ML_e_DL -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_e_DL -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(v))*(1-a1(v)))};
   IL_ML_e_DL -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(v))*(a1(v)))};
   IL_ML_e_DL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-g1(v))*(1))};
; IS_MS_e_DS (Emit: 1,2,4)
   IS_MS_e_DS -> L_L_e_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-g2(v))*((1)*(1)))};
   IS_MS_e_DS -> IS_MS_e_IS  {(g2(v))};
   IS_MS_e_DS -> IS_MS_e_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(v))*(a2(v)))};
   IS_MS_e_DS -> IS_DS_e_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(v))*(a2(v)))};
   IS_MS_e_DS -> IS_MS_e_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_e_DS -> IS_DS_e_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(v))*(1-a2(v)))};
   IS_MS_e_DS -> IS_IS_e_WS  {(b2(t)) * (1-g2(v))};
; IL_ML_IL_WL (Emit: 3)
   IL_ML_IL_WL -> IL_ML_ML_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_ML_IL_WL -> IL_ML_DL_ML  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_ML_IL_WL -> IL_ML_ML_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_ML_IL_WL -> IL_ML_IL_WL  {(b1(u))};
   IL_ML_IL_WL -> IL_DL_ML_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_ML_IL_WL -> IL_DL_DL_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_ML_IL_WL -> IL_DL_ML_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_ML_IL_WL -> IL_ML_DL_DL  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_ML_IL_WL -> IL_DL_DL_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_ML_IL_WL -> IL_IL_WL_WL  {(b1(t)) * (1-b1(u))};
   IL_ML_IL_WL -> IL_Bi[SiL]_WL_WL  {(b1(t)*pS())};
   IL_ML_IL_WL -> IL_ML_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_ML_IL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_ML_IL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_ML_IL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_ML_IL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_ML_IL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_ML_IL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_ML_IL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_ML_IL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_ML_IL_WL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * (1)};
; IL_DL_e_ML (Emit: 1,2,4)
   IL_DL_e_ML -> IL_DL_e_IL  {(b1(v))};
   IL_DL_e_ML -> IL_ML_e_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_DL_e_ML -> IL_DL_e_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_DL_e_ML -> IL_ML_e_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_e_ML -> IL_DL_e_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_e_ML -> IL_IL_e_WL  {(g1(t)) * (1-b1(v))};
   IL_DL_e_ML -> IL_Bi[SiL]_e_WL  {(b1(t)*pS())};
   IL_DL_e_ML -> IL_DL_e_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_e_ML -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_DL_e_ML -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_e_ML -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v)))};
   IL_DL_e_ML -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(a1(v)))};
   IL_DL_e_ML -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-b1(v))*(1))};
; IL_DL_e_DL (Emit: 1,2,4)
   IL_DL_e_DL -> IL_DL_e_IL  {(g1(v))};
   IL_DL_e_DL -> IL_ML_e_ML  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(v))*(a1(v)))};
   IL_DL_e_DL -> IL_DL_e_ML  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(v))*(a1(v)))};
   IL_DL_e_DL -> IL_ML_e_DL  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_e_DL -> IL_DL_e_DL  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_e_DL -> IL_IL_e_WL  {(g1(t)) * (1-g1(v))};
   IL_DL_e_DL -> IL_Bi[SiL]_e_WL  {(b1(t)*pS())};
   IL_DL_e_DL -> IL_DL_e_Bi[SiL]  {(b1(v)*pS())};
   IL_DL_e_DL -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(v))*(a1(v)))};
   IL_DL_e_DL -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_e_DL -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(v))*(1-a1(v)))};
   IL_DL_e_DL -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(v))*(a1(v)))};
   IL_DL_e_DL -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-g1(v))*(1))};
; IS_DS_e_DS (Emit: 1,2,4)
   IS_DS_e_DS -> L_L_e_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-g2(v))*((1)*(1)))};
   IS_DS_e_DS -> IS_DS_e_IS  {(g2(v))};
   IS_DS_e_DS -> IS_MS_e_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(v))*(a2(v)))};
   IS_DS_e_DS -> IS_DS_e_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(v))*(a2(v)))};
   IS_DS_e_DS -> IS_MS_e_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_e_DS -> IS_DS_e_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(v))*(1-a2(v)))};
   IS_DS_e_DS -> IS_IS_e_WS  {(g2(t)) * (1-g2(v))};
; IS_DS_IS_WS (Emit: 3)
   IS_DS_IS_WS -> L_L_L_L  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1))) * ((1)*(1))};
   IS_DS_IS_WS -> IS_MS_MS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * (a2(v))};
   IS_DS_IS_WS -> IS_DS_MS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * (a2(v))};
   IS_DS_IS_WS -> IS_MS_DS_MS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * (a2(v))};
   IS_DS_IS_WS -> IS_MS_MS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u))) * (1-a2(v))};
   IS_DS_IS_WS -> IS_DS_DS_MS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * (a2(v))};
   IS_DS_IS_WS -> IS_DS_MS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u))) * (1-a2(v))};
   IS_DS_IS_WS -> IS_MS_DS_DS  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u))) * (1-a2(v))};
   IS_DS_IS_WS -> IS_DS_DS_DS  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u))) * (1-a2(v))};
   IS_DS_IS_WS -> IS_DS_IS_WS  {(b2(u))};
   IS_DS_IS_WS -> IS_IS_WS_WS  {(g2(t)) * (1-b2(u))};
; IS_DS_IS_e (Emit: 3)
   IS_DS_IS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1)))};
   IS_DS_IS_e -> IS_DS_IS_e  {(b2(u))};
   IS_DS_IS_e -> IS_MS_MS_e  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u)))};
   IS_DS_IS_e -> IS_DS_MS_e  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u)))};
   IS_DS_IS_e -> IS_MS_DS_e  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u)))};
   IS_DS_IS_e -> IS_DS_DS_e  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u)))};
   IS_DS_IS_e -> IS_IS_WS_e  {(g2(t)) * (1-b2(u))};
; IL_DL_IL_e (Emit: 3)
   IL_DL_IL_e -> IL_DL_IL_e  {(b1(u))};
   IL_DL_IL_e -> IL_ML_ML_e  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_DL_IL_e -> IL_DL_ML_e  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_DL_IL_e -> IL_ML_DL_e  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_DL_IL_e -> IL_DL_DL_e  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_DL_IL_e -> IL_IL_WL_e  {(g1(t)) * (1-b1(u))};
   IL_DL_IL_e -> IL_Bi[SiL]_WL_e  {(b1(t)*pS())};
   IL_DL_IL_e -> IL_DL_Bi[SiL]_e  {(b1(u)*pS())};
   IL_DL_IL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_DL_IL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u)))};
   IL_DL_IL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_DL_IL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u)))};
   IL_DL_IL_e -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1))};
; IS_e_e_IS (Emit: 4)
   IS_e_e_IS -> L_e_e_L  {((1-K2())*(1)) * ((1-b2(v))*((1)*(1)))};
   IS_e_e_IS -> IS_e_e_IS  {(b2(v))};
   IS_e_e_IS -> IS_e_e_MS  {(K2()) * ((1-b2(v))*(a2(v)))};
   IS_e_e_IS -> IS_e_e_DS  {(K2()) * ((1-b2(v))*(1-a2(v)))};
; IS_e_MS_MS (Emit: 1,3,4)
   IS_e_MS_MS -> L_e_L_L  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_e_MS_MS -> IS_e_MS_IS  {(b2(v))};
   IS_e_MS_MS -> IS_e_MS_MS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_e_MS_MS -> IS_e_DS_MS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(a2(v)))};
   IS_e_MS_MS -> IS_e_MS_DS  {(K2()) * ((1-b2(u))*(a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_e_MS_MS -> IS_e_DS_DS  {(K2()) * ((1-b2(u))*(1-a2(u))) * ((1-b2(v))*(1-a2(v)))};
   IS_e_MS_MS -> IS_e_IS_WS  {(b2(u)) * (1-b2(v))};
; IS_MS_e_MS (Emit: 1,2,4)
   IS_MS_e_MS -> L_L_e_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-b2(v))*((1)*(1)))};
   IS_MS_e_MS -> IS_MS_e_IS  {(b2(v))};
   IS_MS_e_MS -> IS_MS_e_MS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(v))*(a2(v)))};
   IS_MS_e_MS -> IS_DS_e_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(v))*(a2(v)))};
   IS_MS_e_MS -> IS_MS_e_DS  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_e_MS -> IS_DS_e_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-b2(v))*(1-a2(v)))};
   IS_MS_e_MS -> IS_IS_e_WS  {(b2(t)) * (1-b2(v))};
; IL_e_ML_DL (Emit: 1,3,4)
   IL_e_ML_DL -> IL_e_ML_IL  {(g1(v))};
   IL_e_ML_DL -> IL_e_ML_ML  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_e_ML_DL -> IL_e_DL_ML  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_e_ML_DL -> IL_e_ML_DL  {(K1()) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_e_ML_DL -> IL_e_DL_DL  {(K1()) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_e_ML_DL -> IL_e_IL_WL  {(b1(u)) * (1-g1(v))};
   IL_e_ML_DL -> IL_e_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_e_ML_DL -> IL_e_ML_Bi[SiL]  {(b1(v)*pS())};
   IL_e_ML_DL -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_e_ML_DL -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(a1(v)))};
   IL_e_ML_DL -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_e_ML_DL -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-g1(v))*(1-a1(v)))};
   IL_e_ML_DL -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1)) * ((1-g1(v))*(1))};
; IS_e_DS_DS (Emit: 1,3,4)
   IS_e_DS_DS -> L_e_L_L  {((1-K2())*(1)) * ((1-g2(u))*((1)*(1))) * ((1-g2(v))*((1)*(1)))};
   IS_e_DS_DS -> IS_e_DS_IS  {(g2(v))};
   IS_e_DS_DS -> IS_e_MS_MS  {(K2()) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_e_DS_DS -> IS_e_DS_MS  {(K2()) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(a2(v)))};
   IS_e_DS_DS -> IS_e_MS_DS  {(K2()) * ((1-g2(u))*(a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_e_DS_DS -> IS_e_DS_DS  {(K2()) * ((1-g2(u))*(1-a2(u))) * ((1-g2(v))*(1-a2(v)))};
   IS_e_DS_DS -> IS_e_IS_WS  {(g2(u)) * (1-g2(v))};
; IL_e_e_DL (Emit: 1,4)
   IL_e_e_DL -> IL_e_e_IL  {(g1(v))};
   IL_e_e_DL -> IL_e_e_ML  {(K1()) * ((1-g1(v))*(a1(v)))};
   IL_e_e_DL -> IL_e_e_DL  {(K1()) * ((1-g1(v))*(1-a1(v)))};
   IL_e_e_DL -> IL_e_e_Bi[SiL]  {(b1(v)*pS())};
   IL_e_e_DL -> Bi[SL]_e_e_Bm[SL]  {(K1()*pS()) * ((1-g1(v))*(a1(v)))};
   IL_e_e_DL -> Bi[SL]_e_e_Bp[eL]  {(K1()*pS()) * ((1-g1(v))*(1-a1(v)))};
   IL_e_e_DL -> e_e_e_e  {(1-K1()) * ((1-g1(v))*(1))};
; IL_IL_e_WL (Emit: 2)
   IL_IL_e_WL -> IL_ML_e_ML  {(K1()) * ((1-b1(t))*(a1(t))) * (a1(v))};
   IL_IL_e_WL -> IL_DL_e_ML  {(K1()) * ((1-b1(t))*(1-a1(t))) * (a1(v))};
   IL_IL_e_WL -> IL_ML_e_DL  {(K1()) * ((1-b1(t))*(a1(t))) * (1-a1(v))};
   IL_IL_e_WL -> IL_DL_e_DL  {(K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(v))};
   IL_IL_e_WL -> IL_IL_e_WL  {(b1(t))};
   IL_IL_e_WL -> IL_Bi[SiL]_e_WL  {(b1(t)*pS())};
   IL_IL_e_WL -> Bi[SL]_Bp[eL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(v))};
   IL_IL_e_WL -> Bi[SL]_Bp[eL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(v))};
   IL_IL_e_WL -> Bi[SL]_Bm[SL]_e_Bp[eL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(v))};
   IL_IL_e_WL -> Bi[SL]_Bm[SL]_e_Bm[SL]  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(v))};
   IL_IL_e_WL -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * (1)};
; IL_e_IL_e (Emit: 3)
   IL_e_IL_e -> IL_e_IL_e  {(b1(u))};
   IL_e_IL_e -> IL_e_ML_e  {(K1()) * ((1-b1(u))*(a1(u)))};
   IL_e_IL_e -> IL_e_DL_e  {(K1()) * ((1-b1(u))*(1-a1(u)))};
   IL_e_IL_e -> IL_e_Bi[SiL]_e  {(b1(u)*pS())};
   IL_e_IL_e -> Bi[SL]_e_Bp[eL]_e  {(K1()*pS()) * ((1-b1(u))*(1-a1(u)))};
   IL_e_IL_e -> Bi[SL]_e_Bm[SL]_e  {(K1()*pS()) * ((1-b1(u))*(a1(u)))};
   IL_e_IL_e -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1))};
; IS_e_IS_e (Emit: 3)
   IS_e_IS_e -> L_e_L_e  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1)))};
   IS_e_IS_e -> IS_e_IS_e  {(b2(u))};
   IS_e_IS_e -> IS_e_MS_e  {(K2()) * ((1-b2(u))*(a2(u)))};
   IS_e_IS_e -> IS_e_DS_e  {(K2()) * ((1-b2(u))*(1-a2(u)))};
; IS_DS_MS_e (Emit: 1,2,3)
   IS_DS_MS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-b2(u))*((1)*(1)))};
   IS_DS_MS_e -> IS_DS_IS_e  {(b2(u))};
   IS_DS_MS_e -> IS_MS_MS_e  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(a2(u)))};
   IS_DS_MS_e -> IS_DS_MS_e  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(a2(u)))};
   IS_DS_MS_e -> IS_MS_DS_e  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-b2(u))*(1-a2(u)))};
   IS_DS_MS_e -> IS_DS_DS_e  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-b2(u))*(1-a2(u)))};
   IS_DS_MS_e -> IS_IS_WS_e  {(g2(t)) * (1-b2(u))};
; IL_ML_DL_e (Emit: 1,2,3)
   IL_ML_DL_e -> IL_ML_IL_e  {(g1(u))};
   IL_ML_DL_e -> IL_ML_ML_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u)))};
   IL_ML_DL_e -> IL_DL_ML_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u)))};
   IL_ML_DL_e -> IL_ML_DL_e  {(K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u)))};
   IL_ML_DL_e -> IL_DL_DL_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u)))};
   IL_ML_DL_e -> IL_IL_WL_e  {(b1(t)) * (1-g1(u))};
   IL_ML_DL_e -> IL_Bi[SiL]_WL_e  {(b1(t)*pS())};
   IL_ML_DL_e -> IL_ML_Bi[SiL]_e  {(b1(u)*pS())};
   IL_ML_DL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u)))};
   IL_ML_DL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u)))};
   IL_ML_DL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u)))};
   IL_ML_DL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u)))};
   IL_ML_DL_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * ((1-g1(u))*(1))};
; IL_DL_DL_e (Emit: 1,2,3)
   IL_DL_DL_e -> IL_DL_IL_e  {(g1(u))};
   IL_DL_DL_e -> IL_ML_ML_e  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u)))};
   IL_DL_DL_e -> IL_DL_ML_e  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u)))};
   IL_DL_DL_e -> IL_ML_DL_e  {(K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u)))};
   IL_DL_DL_e -> IL_DL_DL_e  {(K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u)))};
   IL_DL_DL_e -> IL_IL_WL_e  {(g1(t)) * (1-g1(u))};
   IL_DL_DL_e -> IL_Bi[SiL]_WL_e  {(b1(t)*pS())};
   IL_DL_DL_e -> IL_DL_Bi[SiL]_e  {(b1(u)*pS())};
   IL_DL_DL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u)))};
   IL_DL_DL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u)))};
   IL_DL_DL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u)))};
   IL_DL_DL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u)))};
   IL_DL_DL_e -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1)) * ((1-g1(u))*(1))};
; IL_IL_WL_e (Emit: 2)
   IL_IL_WL_e -> IL_ML_ML_e  {(K1()) * ((1-b1(t))*(a1(t))) * (a1(u))};
   IL_IL_WL_e -> IL_DL_ML_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u))};
   IL_IL_WL_e -> IL_ML_DL_e  {(K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u))};
   IL_IL_WL_e -> IL_DL_DL_e  {(K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u))};
   IL_IL_WL_e -> IL_IL_WL_e  {(b1(t))};
   IL_IL_WL_e -> IL_Bi[SiL]_WL_e  {(b1(t)*pS())};
   IL_IL_WL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u))};
   IL_IL_WL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u))};
   IL_IL_WL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u))};
   IL_IL_WL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e  {(K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u))};
   IL_IL_WL_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1)) * (1)};
; IL_e_e_IL (Emit: 4)
   IL_e_e_IL -> IL_e_e_IL  {(b1(v))};
   IL_e_e_IL -> IL_e_e_ML  {(K1()) * ((1-b1(v))*(a1(v)))};
   IL_e_e_IL -> IL_e_e_DL  {(K1()) * ((1-b1(v))*(1-a1(v)))};
   IL_e_e_IL -> IL_e_e_Bi[SiL]  {(b1(v)*pS())};
   IL_e_e_IL -> Bi[SL]_e_e_Bm[SL]  {(K1()*pS()) * ((1-b1(v))*(a1(v)))};
   IL_e_e_IL -> Bi[SL]_e_e_Bp[eL]  {(K1()*pS()) * ((1-b1(v))*(1-a1(v)))};
   IL_e_e_IL -> e_e_e_e  {(1-K1()) * ((1-b1(v))*(1))};
; IS_e_e_DS (Emit: 1,4)
   IS_e_e_DS -> L_e_e_L  {((1-K2())*(1)) * ((1-g2(v))*((1)*(1)))};
   IS_e_e_DS -> IS_e_e_IS  {(g2(v))};
   IS_e_e_DS -> IS_e_e_MS  {(K2()) * ((1-g2(v))*(a2(v)))};
   IS_e_e_DS -> IS_e_e_DS  {(K2()) * ((1-g2(v))*(1-a2(v)))};
; IS_e_IS_WS (Emit: 3)
   IS_e_IS_WS -> L_e_L_L  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1))) * ((1)*(1))};
   IS_e_IS_WS -> IS_e_MS_MS  {(K2()) * ((1-b2(u))*(a2(u))) * (a2(v))};
   IS_e_IS_WS -> IS_e_DS_MS  {(K2()) * ((1-b2(u))*(1-a2(u))) * (a2(v))};
   IS_e_IS_WS -> IS_e_MS_DS  {(K2()) * ((1-b2(u))*(a2(u))) * (1-a2(v))};
   IS_e_IS_WS -> IS_e_DS_DS  {(K2()) * ((1-b2(u))*(1-a2(u))) * (1-a2(v))};
   IS_e_IS_WS -> IS_e_IS_WS  {(b2(u))};
; IL_e_IL_WL (Emit: 3)
   IL_e_IL_WL -> IL_e_ML_ML  {(K1()) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_e_IL_WL -> IL_e_DL_ML  {(K1()) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_e_IL_WL -> IL_e_ML_DL  {(K1()) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_e_IL_WL -> IL_e_DL_DL  {(K1()) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_e_IL_WL -> IL_e_IL_WL  {(b1(u))};
   IL_e_IL_WL -> IL_e_Bi[SiL]_WL  {(b1(u)*pS())};
   IL_e_IL_WL -> Bi[SL]_e_Bp[eL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * (a1(v))};
   IL_e_IL_WL -> Bi[SL]_e_Bm[SL]_Bm[SL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * (a1(v))};
   IL_e_IL_WL -> Bi[SL]_e_Bp[eL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(1-a1(u))) * (1-a1(v))};
   IL_e_IL_WL -> Bi[SL]_e_Bm[SL]_Bp[eL]  {(K1()*pS()) * ((1-b1(u))*(a1(u))) * (1-a1(v))};
   IL_e_IL_WL -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1)) * (1)};
; IS_IS_WS_WS (Emit: 2)
   IS_IS_WS_WS -> L_L_L_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1)*(1)) * ((1)*(1))};
   IS_IS_WS_WS -> IS_MS_MS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * (a2(u)) * (a2(v))};
   IS_IS_WS_WS -> IS_DS_MS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (a2(u)) * (a2(v))};
   IS_IS_WS_WS -> IS_MS_DS_MS  {(K2()) * ((1-b2(t))*(a2(t))) * (1-a2(u)) * (a2(v))};
   IS_IS_WS_WS -> IS_MS_MS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * (a2(u)) * (1-a2(v))};
   IS_IS_WS_WS -> IS_DS_DS_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (1-a2(u)) * (a2(v))};
   IS_IS_WS_WS -> IS_DS_MS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (a2(u)) * (1-a2(v))};
   IS_IS_WS_WS -> IS_MS_DS_DS  {(K2()) * ((1-b2(t))*(a2(t))) * (1-a2(u)) * (1-a2(v))};
   IS_IS_WS_WS -> IS_DS_DS_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (1-a2(u)) * (1-a2(v))};
   IS_IS_WS_WS -> IS_IS_WS_WS  {(b2(t))};
; IS_IS_e_WS (Emit: 2)
   IS_IS_e_WS -> L_L_e_L  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1)*(1))};
   IS_IS_e_WS -> IS_MS_e_MS  {(K2()) * ((1-b2(t))*(a2(t))) * (a2(v))};
   IS_IS_e_WS -> IS_DS_e_MS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (a2(v))};
   IS_IS_e_WS -> IS_MS_e_DS  {(K2()) * ((1-b2(t))*(a2(t))) * (1-a2(v))};
   IS_IS_e_WS -> IS_DS_e_DS  {(K2()) * ((1-b2(t))*(1-a2(t))) * (1-a2(v))};
   IS_IS_e_WS -> IS_IS_e_WS  {(b2(t))};
; IS_e_MS_e (Emit: 1,3)
   IS_e_MS_e -> L_e_L_e  {((1-K2())*(1)) * ((1-b2(u))*((1)*(1)))};
   IS_e_MS_e -> IS_e_IS_e  {(b2(u))};
   IS_e_MS_e -> IS_e_MS_e  {(K2()) * ((1-b2(u))*(a2(u)))};
   IS_e_MS_e -> IS_e_DS_e  {(K2()) * ((1-b2(u))*(1-a2(u)))};
; IL_e_ML_e (Emit: 1,3)
   IL_e_ML_e -> IL_e_IL_e  {(b1(u))};
   IL_e_ML_e -> IL_e_ML_e  {(K1()) * ((1-b1(u))*(a1(u)))};
   IL_e_ML_e -> IL_e_DL_e  {(K1()) * ((1-b1(u))*(1-a1(u)))};
   IL_e_ML_e -> IL_e_Bi[SiL]_e  {(b1(u)*pS())};
   IL_e_ML_e -> Bi[SL]_e_Bp[eL]_e  {(K1()*pS()) * ((1-b1(u))*(1-a1(u)))};
   IL_e_ML_e -> Bi[SL]_e_Bm[SL]_e  {(K1()*pS()) * ((1-b1(u))*(a1(u)))};
   IL_e_ML_e -> e_e_e_e  {(1-K1()) * ((1-b1(u))*(1))};
; IS_MS_DS_e (Emit: 1,2,3)
   IS_MS_DS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1-g2(u))*((1)*(1)))};
   IS_MS_DS_e -> IS_MS_IS_e  {(g2(u))};
   IS_MS_DS_e -> IS_MS_MS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(a2(u)))};
   IS_MS_DS_e -> IS_DS_MS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u)))};
   IS_MS_DS_e -> IS_MS_DS_e  {(K2()) * ((1-b2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u)))};
   IS_MS_DS_e -> IS_DS_DS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u)))};
   IS_MS_DS_e -> IS_IS_WS_e  {(b2(t)) * (1-g2(u))};
; IS_DS_DS_e (Emit: 1,2,3)
   IS_DS_DS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1))) * ((1-g2(u))*((1)*(1)))};
   IS_DS_DS_e -> IS_DS_IS_e  {(g2(u))};
   IS_DS_DS_e -> IS_MS_MS_e  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(a2(u)))};
   IS_DS_DS_e -> IS_DS_MS_e  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(a2(u)))};
   IS_DS_DS_e -> IS_MS_DS_e  {(K2()) * ((1-g2(t))*(a2(t))) * ((1-g2(u))*(1-a2(u)))};
   IS_DS_DS_e -> IS_DS_DS_e  {(K2()) * ((1-g2(t))*(1-a2(t))) * ((1-g2(u))*(1-a2(u)))};
   IS_DS_DS_e -> IS_IS_WS_e  {(g2(t)) * (1-g2(u))};
; IS_e_DS_e (Emit: 1,3)
   IS_e_DS_e -> L_e_L_e  {((1-K2())*(1)) * ((1-g2(u))*((1)*(1)))};
   IS_e_DS_e -> IS_e_IS_e  {(g2(u))};
   IS_e_DS_e -> IS_e_MS_e  {(K2()) * ((1-g2(u))*(a2(u)))};
   IS_e_DS_e -> IS_e_DS_e  {(K2()) * ((1-g2(u))*(1-a2(u)))};
; IL_e_DL_e (Emit: 1,3)
   IL_e_DL_e -> IL_e_IL_e  {(g1(u))};
   IL_e_DL_e -> IL_e_ML_e  {(K1()) * ((1-g1(u))*(a1(u)))};
   IL_e_DL_e -> IL_e_DL_e  {(K1()) * ((1-g1(u))*(1-a1(u)))};
   IL_e_DL_e -> IL_e_Bi[SiL]_e  {(b1(u)*pS())};
   IL_e_DL_e -> Bi[SL]_e_Bp[eL]_e  {(K1()*pS()) * ((1-g1(u))*(1-a1(u)))};
   IL_e_DL_e -> Bi[SL]_e_Bm[SL]_e  {(K1()*pS()) * ((1-g1(u))*(a1(u)))};
   IL_e_DL_e -> e_e_e_e  {(1-K1()) * ((1-g1(u))*(1))};
; IS_IS_WS_e (Emit: 2)
   IS_IS_WS_e -> L_L_L_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1))) * ((1)*(1))};
   IS_IS_WS_e -> IS_MS_MS_e  {(K2()) * ((1-b2(t))*(a2(t))) * (a2(u))};
   IS_IS_WS_e -> IS_DS_MS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * (a2(u))};
   IS_IS_WS_e -> IS_MS_DS_e  {(K2()) * ((1-b2(t))*(a2(t))) * (1-a2(u))};
   IS_IS_WS_e -> IS_DS_DS_e  {(K2()) * ((1-b2(t))*(1-a2(t))) * (1-a2(u))};
   IS_IS_WS_e -> IS_IS_WS_e  {(b2(t))};
; IS_IS_e_e (Emit: 2)
   IS_IS_e_e -> L_L_e_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1)))};
   IS_IS_e_e -> IS_IS_e_e  {(b2(t))};
   IS_IS_e_e -> IS_MS_e_e  {(K2()) * ((1-b2(t))*(a2(t)))};
   IS_IS_e_e -> IS_DS_e_e  {(K2()) * ((1-b2(t))*(1-a2(t)))};
; IL_IL_e_e (Emit: 2)
   IL_IL_e_e -> IL_IL_e_e  {(b1(t))};
   IL_IL_e_e -> IL_ML_e_e  {(K1()) * ((1-b1(t))*(a1(t)))};
   IL_IL_e_e -> IL_DL_e_e  {(K1()) * ((1-b1(t))*(1-a1(t)))};
   IL_IL_e_e -> IL_Bi[SiL]_e_e  {(b1(t)*pS())};
   IL_IL_e_e -> Bi[SL]_Bp[eL]_e_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t)))};
   IL_IL_e_e -> Bi[SL]_Bm[SL]_e_e  {(K1()*pS()) * ((1-b1(t))*(a1(t)))};
   IL_IL_e_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1))};
; IS_MS_e_e (Emit: 1,2)
   IS_MS_e_e -> L_L_e_e  {((1-K2())*(1)) * ((1-b2(t))*((1)*(1)))};
   IS_MS_e_e -> IS_IS_e_e  {(b2(t))};
   IS_MS_e_e -> IS_MS_e_e  {(K2()) * ((1-b2(t))*(a2(t)))};
   IS_MS_e_e -> IS_DS_e_e  {(K2()) * ((1-b2(t))*(1-a2(t)))};
; IL_ML_e_e (Emit: 1,2)
   IL_ML_e_e -> IL_IL_e_e  {(b1(t))};
   IL_ML_e_e -> IL_ML_e_e  {(K1()) * ((1-b1(t))*(a1(t)))};
   IL_ML_e_e -> IL_DL_e_e  {(K1()) * ((1-b1(t))*(1-a1(t)))};
   IL_ML_e_e -> IL_Bi[SiL]_e_e  {(b1(t)*pS())};
   IL_ML_e_e -> Bi[SL]_Bp[eL]_e_e  {(K1()*pS()) * ((1-b1(t))*(1-a1(t)))};
   IL_ML_e_e -> Bi[SL]_Bm[SL]_e_e  {(K1()*pS()) * ((1-b1(t))*(a1(t)))};
   IL_ML_e_e -> e_e_e_e  {(1-K1()) * ((1-b1(t))*(1))};
; IL_DL_e_e (Emit: 1,2)
   IL_DL_e_e -> IL_IL_e_e  {(g1(t))};
   IL_DL_e_e -> IL_ML_e_e  {(K1()) * ((1-g1(t))*(a1(t)))};
   IL_DL_e_e -> IL_DL_e_e  {(K1()) * ((1-g1(t))*(1-a1(t)))};
   IL_DL_e_e -> IL_Bi[SiL]_e_e  {(b1(t)*pS())};
   IL_DL_e_e -> Bi[SL]_Bp[eL]_e_e  {(K1()*pS()) * ((1-g1(t))*(1-a1(t)))};
   IL_DL_e_e -> Bi[SL]_Bm[SL]_e_e  {(K1()*pS()) * ((1-g1(t))*(a1(t)))};
   IL_DL_e_e -> e_e_e_e  {(1-K1()) * ((1-g1(t))*(1))};
; IS_DS_e_e (Emit: 1,2)
   IS_DS_e_e -> L_L_e_e  {((1-K2())*(1)) * ((1-g2(t))*((1)*(1)))};
   IS_DS_e_e -> IS_IS_e_e  {(g2(t))};
   IS_DS_e_e -> IS_MS_e_e  {(K2()) * ((1-g2(t))*(a2(t)))};
   IS_DS_e_e -> IS_DS_e_e  {(K2()) * ((1-g2(t))*(1-a2(t)))};
; IL_e_e_e (Emit: 1)
   IL_e_e_e -> IL_e_e_e  {(K1())};
   IL_e_e_e -> Bi[SL]_e_e_e  {(K1()*pS())};
   IL_e_e_e -> e_e_e_e  {(1-K1())};
; IS_e_e_e (Emit: 1)
   IS_e_e_e -> L_e_e_e  {((1-K2())*(1))};
   IS_e_e_e -> IS_e_e_e  {(K2())};
; e_e_e_Li (Start)
   e_e_e_Li -> e_e_e_ILi  {(K1())};
   e_e_e_Li -> e_e_e_Bi[SiLi]  {(K1()*pS())};
   e_e_e_Li -> e_e_e_e  {(1-K1())};
; e_e_e_Si (Start)
   e_e_e_Si -> e_e_e_Li  {((1-K2())*(1))};
   e_e_e_Si -> e_e_e_ISi  {(K2())};
; e_e_e_ILi (Emit: 4)
   e_e_e_ILi -> e_e_e_ILi  {(K1())};
   e_e_e_ILi -> e_e_e_Bi[SiLi]  {(K1()*pS())};
   e_e_e_ILi -> e_e_e_e  {(1-K1())};
; e_e_e_ISi (Emit: 4)
   e_e_e_ISi -> e_e_e_Li  {((1-K2())*(1))};
   e_e_e_ISi -> e_e_e_ISi  {(K2())};
; e_e_Li_e (Start)
   e_e_Li_e -> e_e_ILi_e  {(K1())};
   e_e_Li_e -> e_e_Bi[SiLi]_e  {(K1()*pS())};
   e_e_Li_e -> e_e_e_e  {(1-K1())};
; e_e_Si_e (Start)
   e_e_Si_e -> e_e_Li_e  {((1-K2())*(1))};
   e_e_Si_e -> e_e_ISi_e  {(K2())};
; e_e_ILi_e (Emit: 3)
   e_e_ILi_e -> e_e_ILi_e  {(K1())};
   e_e_ILi_e -> e_e_Bi[SiLi]_e  {(K1()*pS())};
   e_e_ILi_e -> e_e_e_e  {(1-K1())};
; e_e_ISi_e (Emit: 3)
   e_e_ISi_e -> e_e_Li_e  {((1-K2())*(1))};
   e_e_ISi_e -> e_e_ISi_e  {(K2())};
; e_Si_e_e (Start)
   e_Si_e_e -> e_Li_e_e  {((1-K2())*(1))};
   e_Si_e_e -> e_ISi_e_e  {(K2())};
; e_Li_e_e (Start)
   e_Li_e_e -> e_ILi_e_e  {(K1())};
   e_Li_e_e -> e_Bi[SiLi]_e_e  {(K1()*pS())};
   e_Li_e_e -> e_e_e_e  {(1-K1())};
; e_ISi_e_e (Emit: 2)
   e_ISi_e_e -> e_Li_e_e  {((1-K2())*(1))};
   e_ISi_e_e -> e_ISi_e_e  {(K2())};
; e_ILi_e_e (Emit: 2)
   e_ILi_e_e -> e_ILi_e_e  {(K1())};
   e_ILi_e_e -> e_Bi[SiLi]_e_e  {(K1()*pS())};
   e_ILi_e_e -> e_e_e_e  {(1-K1())};

;; Bifurcations
; L_L_Bi[SiL]_WL (Bifurc)
;  (left):
      L_L_Bi[SiL]_WL -> e_e_Si_e {(1)};
;  (center):
      L_L_Bi[SiL]_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * (1))};
      L_L_Bi[SiL]_WL -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      L_L_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      L_L_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      L_L_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      L_L_Bi[SiL]_WL -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      L_L_Bi[SiL]_WL -> L_L_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      L_L_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      L_L_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      L_L_Bi[SiL]_WL -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      L_L_Bi[SiL]_WL -> L_IL_WL_WL {((1)) * ((b1(t)) * (1-b1(u)))};
      L_L_Bi[SiL]_WL -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      L_L_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      L_L_Bi[SiL]_WL -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      L_L_Bi[SiL]_WL -> L_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      L_L_Bi[SiL]_WL -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      L_L_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      L_L_Bi[SiL]_WL -> L_L_IL_WL {((1)) * ((b1(u)))};
      L_L_Bi[SiL]_WL -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      L_L_Bi[SiL]_WL -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      L_L_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};

; L_L_L_Bi[SiL] (Bifurc)
;  (left):
      L_L_L_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      L_L_L_Bi[SiL] -> L_L_L_L {(1)};

; L_Bi[SiL]_WL_WL (Bifurc)
;  (left):
      L_Bi[SiL]_WL_WL -> e_Si_e_e {(1)};
;  (center):
      L_Bi[SiL]_WL_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * (1) * (1))};
      L_Bi[SiL]_WL_WL -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (1-a1(v)))};
      L_Bi[SiL]_WL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (1-a1(v)))};
      L_Bi[SiL]_WL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (1-a1(v)))};
      L_Bi[SiL]_WL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v)))};
      L_Bi[SiL]_WL_WL -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (1-a1(v)))};
      L_Bi[SiL]_WL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (a1(v)))};
      L_Bi[SiL]_WL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (a1(v)))};
      L_Bi[SiL]_WL_WL -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (a1(v)))};
      L_Bi[SiL]_WL_WL -> L_IL_WL_WL {((1)) * ((b1(t)))};
      L_Bi[SiL]_WL_WL -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (a1(v)))};
      L_Bi[SiL]_WL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (1-a1(v)))};
      L_Bi[SiL]_WL_WL -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (1-a1(v)))};
      L_Bi[SiL]_WL_WL -> L_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      L_Bi[SiL]_WL_WL -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (1-a1(v)))};
      L_Bi[SiL]_WL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (1-a1(v)))};
      L_Bi[SiL]_WL_WL -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (a1(v)))};
      L_Bi[SiL]_WL_WL -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v)))};
      L_Bi[SiL]_WL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (a1(v)))};

; L_e_Bi[SiL]_WL (Bifurc)
;  (left):
      L_e_Bi[SiL]_WL -> e_e_Si_e {(1)};
;  (center):
      L_e_Bi[SiL]_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(u))*(1)) * (1))};
      L_e_Bi[SiL]_WL -> Bi[SL]_e_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      L_e_Bi[SiL]_WL -> L_e_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      L_e_Bi[SiL]_WL -> L_e_IL_WL {((1)) * ((b1(u)))};
      L_e_Bi[SiL]_WL -> IL_e_ML_DL {((1)) * ((K1()) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      L_e_Bi[SiL]_WL -> IL_e_ML_ML {((1)) * ((K1()) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      L_e_Bi[SiL]_WL -> Bi[SL]_e_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      L_e_Bi[SiL]_WL -> IL_e_DL_DL {((1)) * ((K1()) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      L_e_Bi[SiL]_WL -> IL_e_DL_ML {((1)) * ((K1()) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      L_e_Bi[SiL]_WL -> Bi[SL]_e_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      L_e_Bi[SiL]_WL -> Bi[SL]_e_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};

; L_e_L_Bi[SiL] (Bifurc)
;  (left):
      L_e_L_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      L_e_L_Bi[SiL] -> L_e_L_L {(1)};

; L_L_e_Bi[SiL] (Bifurc)
;  (left):
      L_L_e_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      L_L_e_Bi[SiL] -> L_L_e_L {(1)};

; L_Bi[SiL]_WL_e (Bifurc)
;  (left):
      L_Bi[SiL]_WL_e -> e_Si_e_e {(1)};
;  (center):
      L_Bi[SiL]_WL_e -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * (1))};
      L_Bi[SiL]_WL_e -> IL_ML_ML_e {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (a1(u)))};
      L_Bi[SiL]_WL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)))};
      L_Bi[SiL]_WL_e -> L_IL_WL_e {((1)) * ((b1(t)))};
      L_Bi[SiL]_WL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)))};
      L_Bi[SiL]_WL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)))};
      L_Bi[SiL]_WL_e -> IL_DL_DL_e {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)))};
      L_Bi[SiL]_WL_e -> L_Bi[SiL]_WL_e {((1)) * ((b1(t)*pS()))};
      L_Bi[SiL]_WL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)))};
      L_Bi[SiL]_WL_e -> IL_DL_ML_e {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)))};
      L_Bi[SiL]_WL_e -> IL_ML_DL_e {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)))};

; L_L_Bi[SiL]_e (Bifurc)
;  (left):
      L_L_Bi[SiL]_e -> e_e_Si_e {(1)};
;  (center):
      L_L_Bi[SiL]_e -> L_L_L_e {(1)};

; L_Bi[SiL]_e_WL (Bifurc)
;  (left):
      L_Bi[SiL]_e_WL -> e_Si_e_e {(1)};
;  (center):
      L_Bi[SiL]_e_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * (1))};
      L_Bi[SiL]_e_WL -> Bi[SL]_Bp[eL]_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(v)))};
      L_Bi[SiL]_e_WL -> Bi[SL]_Bp[eL]_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(v)))};
      L_Bi[SiL]_e_WL -> IL_DL_e_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (a1(v)))};
      L_Bi[SiL]_e_WL -> L_Bi[SiL]_e_WL {((1)) * ((b1(t)*pS()))};
      L_Bi[SiL]_e_WL -> IL_ML_e_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (a1(v)))};
      L_Bi[SiL]_e_WL -> Bi[SL]_Bm[SL]_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(v)))};
      L_Bi[SiL]_e_WL -> Bi[SL]_Bm[SL]_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(v)))};
      L_Bi[SiL]_e_WL -> IL_ML_e_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (1-a1(v)))};
      L_Bi[SiL]_e_WL -> L_IL_e_WL {((1)) * ((b1(t)))};
      L_Bi[SiL]_e_WL -> IL_DL_e_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(v)))};

; L_e_e_Bi[SiL] (Bifurc)
;  (left):
      L_e_e_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      L_e_e_Bi[SiL] -> L_e_e_L {(1)};

; L_e_Bi[SiL]_e (Bifurc)
;  (left):
      L_e_Bi[SiL]_e -> e_e_Si_e {(1)};
;  (center):
      L_e_Bi[SiL]_e -> L_e_L_e {(1)};

; L_Bi[SiL]_e_e (Bifurc)
;  (left):
      L_Bi[SiL]_e_e -> e_Si_e_e {(1)};
;  (center):
      L_Bi[SiL]_e_e -> L_L_e_e {(1)};

; IL_ML_Bi[SiL]_WL (Bifurc)
;  (left):
      IL_ML_Bi[SiL]_WL -> e_e_Si_e {(1)};
;  (center):
      IL_ML_Bi[SiL]_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * (1))};
      IL_ML_Bi[SiL]_WL -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_ML_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_ML_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_ML_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_ML_Bi[SiL]_WL -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_ML_Bi[SiL]_WL -> IL_ML_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_ML_Bi[SiL]_WL -> IL_ML_IL_WL {((1)) * ((b1(u)))};
      IL_ML_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_ML_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      IL_ML_Bi[SiL]_WL -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_ML_Bi[SiL]_WL -> IL_IL_WL_WL {((1)) * ((b1(t)) * (1-b1(u)))};
      IL_ML_Bi[SiL]_WL -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      IL_ML_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_ML_Bi[SiL]_WL -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_ML_Bi[SiL]_WL -> IL_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      IL_ML_Bi[SiL]_WL -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_ML_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_ML_Bi[SiL]_WL -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      IL_ML_Bi[SiL]_WL -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_ML_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};

; IL_ML_ML_Bi[SiL] (Bifurc)
;  (left):
      IL_ML_ML_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_ML_ML_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1)))};
      IL_ML_ML_Bi[SiL] -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_ML_Bi[SiL] -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_ML_Bi[SiL] -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_ML_ML_Bi[SiL] {((1)) * ((b1(v)*pS()))};
      IL_ML_ML_Bi[SiL] -> IL_ML_ML_IL {((1)) * ((b1(v)))};
      IL_ML_ML_Bi[SiL] -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_ML_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_ML_ML_Bi[SiL] -> IL_ML_IL_WL {((1)) * ((b1(u)) * (1-b1(v)))};
      IL_ML_ML_Bi[SiL] -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_ML_Bi[SiL] -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_IL_WL_WL {((1)) * ((b1(t)) * (1-b1(u)) * (1-b1(v)))};
      IL_ML_ML_Bi[SiL] -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_ML_Bi[SiL] -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      IL_ML_ML_Bi[SiL] -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_ML_Bi[SiL] -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_ML_Bi[SiL] -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_ML_Bi[SiL] -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};

; IL_DL_ML_Bi[SiL] (Bifurc)
;  (left):
      IL_DL_ML_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_DL_ML_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1)) * ((1-b1(v))*(1)))};
      IL_DL_ML_Bi[SiL] -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_ML_Bi[SiL] -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_ML_Bi[SiL] -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_ML_Bi[SiL] -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_DL_ML_Bi[SiL] {((1)) * ((b1(v)*pS()))};
      IL_DL_ML_Bi[SiL] -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_ML_Bi[SiL] -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_DL_ML_IL {((1)) * ((b1(v)))};
      IL_DL_ML_Bi[SiL] -> IL_IL_WL_WL {((1)) * ((g1(t)) * (1-b1(u)) * (1-b1(v)))};
      IL_DL_ML_Bi[SiL] -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_ML_Bi[SiL] -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      IL_DL_ML_Bi[SiL] -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_ML_Bi[SiL] -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_DL_IL_WL {((1)) * ((b1(u)) * (1-b1(v)))};
      IL_DL_ML_Bi[SiL] -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_ML_Bi[SiL] -> IL_DL_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_DL_ML_Bi[SiL] -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};

; IL_ML_DL_Bi[SiL] (Bifurc)
;  (left):
      IL_ML_DL_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_ML_DL_Bi[SiL] -> IL_ML_DL_IL {((1)) * ((b1(v)))};
      IL_ML_DL_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * ((1-g1(u))*(1)) * ((1-b1(v))*(1)))};
      IL_ML_DL_Bi[SiL] -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_DL_Bi[SiL] -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_DL_Bi[SiL] -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_DL_Bi[SiL] -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_DL_Bi[SiL] -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_DL_Bi[SiL] -> IL_ML_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_ML_DL_Bi[SiL] -> IL_ML_IL_WL {((1)) * ((g1(u)) * (1-b1(v)))};
      IL_ML_DL_Bi[SiL] -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_DL_Bi[SiL] -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_DL_Bi[SiL] -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_DL_Bi[SiL] -> IL_IL_WL_WL {((1)) * ((b1(t)) * (1-g1(u)) * (1-b1(v)))};
      IL_ML_DL_Bi[SiL] -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_DL_Bi[SiL] -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_DL_Bi[SiL] -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_DL_Bi[SiL] -> IL_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      IL_ML_DL_Bi[SiL] -> IL_ML_DL_Bi[SiL] {((1)) * ((b1(v)*pS()))};
      IL_ML_DL_Bi[SiL] -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_DL_Bi[SiL] -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_DL_Bi[SiL] -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_DL_Bi[SiL] -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_ML_DL_Bi[SiL] -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};

; IL_DL_Bi[SiL]_WL (Bifurc)
;  (left):
      IL_DL_Bi[SiL]_WL -> e_e_Si_e {(1)};
;  (center):
      IL_DL_Bi[SiL]_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1)) * (1))};
      IL_DL_Bi[SiL]_WL -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_DL_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_DL_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_DL_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_DL_Bi[SiL]_WL -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_DL_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_DL_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      IL_DL_Bi[SiL]_WL -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_DL_Bi[SiL]_WL -> IL_IL_WL_WL {((1)) * ((g1(t)) * (1-b1(u)))};
      IL_DL_Bi[SiL]_WL -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      IL_DL_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_DL_Bi[SiL]_WL -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_DL_Bi[SiL]_WL -> IL_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      IL_DL_Bi[SiL]_WL -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_DL_Bi[SiL]_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_DL_Bi[SiL]_WL -> IL_DL_IL_WL {((1)) * ((b1(u)))};
      IL_DL_Bi[SiL]_WL -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      IL_DL_Bi[SiL]_WL -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_DL_Bi[SiL]_WL -> IL_DL_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_DL_Bi[SiL]_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};

; IL_DL_DL_Bi[SiL] (Bifurc)
;  (left):
      IL_DL_DL_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_DL_DL_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-g1(t))*(1)) * ((1-g1(u))*(1)) * ((1-b1(v))*(1)))};
      IL_DL_DL_Bi[SiL] -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_DL_Bi[SiL] -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_DL_Bi[SiL] -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_DL_Bi[SiL] -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_DL_Bi[SiL] -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_DL_Bi[SiL] -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_DL_Bi[SiL] -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_DL_Bi[SiL] -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_DL_Bi[SiL] -> IL_IL_WL_WL {((1)) * ((g1(t)) * (1-g1(u)) * (1-b1(v)))};
      IL_DL_DL_Bi[SiL] -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_DL_Bi[SiL] -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_DL_Bi[SiL] -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_DL_Bi[SiL] -> IL_DL_DL_IL {((1)) * ((b1(v)))};
      IL_DL_DL_Bi[SiL] -> IL_DL_DL_Bi[SiL] {((1)) * ((b1(v)*pS()))};
      IL_DL_DL_Bi[SiL] -> IL_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      IL_DL_DL_Bi[SiL] -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_DL_Bi[SiL] -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_DL_Bi[SiL] -> IL_DL_IL_WL {((1)) * ((g1(u)) * (1-b1(v)))};
      IL_DL_DL_Bi[SiL] -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_DL_Bi[SiL] -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_DL_DL_Bi[SiL] -> IL_DL_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_DL_DL_Bi[SiL] -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};

; IL_Bi[SiL]_WL_WL (Bifurc)
;  (left):
      IL_Bi[SiL]_WL_WL -> e_Si_e_e {(1)};
;  (center):
      IL_Bi[SiL]_WL_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * (1) * (1))};
      IL_Bi[SiL]_WL_WL -> IL_DL_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (1-a1(v)))};
      IL_Bi[SiL]_WL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (1-a1(v)))};
      IL_Bi[SiL]_WL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (1-a1(v)))};
      IL_Bi[SiL]_WL_WL -> Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v)))};
      IL_Bi[SiL]_WL_WL -> IL_DL_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (1-a1(v)))};
      IL_Bi[SiL]_WL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (a1(v)))};
      IL_Bi[SiL]_WL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (a1(v)))};
      IL_Bi[SiL]_WL_WL -> IL_DL_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (a1(v)))};
      IL_Bi[SiL]_WL_WL -> IL_IL_WL_WL {((1)) * ((b1(t)))};
      IL_Bi[SiL]_WL_WL -> IL_DL_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (a1(v)))};
      IL_Bi[SiL]_WL_WL -> Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)) * (1-a1(v)))};
      IL_Bi[SiL]_WL_WL -> IL_ML_ML_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (1-a1(v)))};
      IL_Bi[SiL]_WL_WL -> IL_Bi[SiL]_WL_WL {((1)) * ((b1(t)*pS()))};
      IL_Bi[SiL]_WL_WL -> IL_ML_DL_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (1-a1(v)))};
      IL_Bi[SiL]_WL_WL -> Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (1-a1(v)))};
      IL_Bi[SiL]_WL_WL -> IL_ML_DL_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)) * (a1(v)))};
      IL_Bi[SiL]_WL_WL -> IL_ML_ML_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (a1(u)) * (a1(v)))};
      IL_Bi[SiL]_WL_WL -> Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)) * (a1(v)))};

; IL_e_Bi[SiL]_WL (Bifurc)
;  (left):
      IL_e_Bi[SiL]_WL -> e_e_Si_e {(1)};
;  (center):
      IL_e_Bi[SiL]_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(u))*(1)) * (1))};
      IL_e_Bi[SiL]_WL -> Bi[SL]_e_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      IL_e_Bi[SiL]_WL -> IL_e_IL_WL {((1)) * ((b1(u)))};
      IL_e_Bi[SiL]_WL -> IL_e_ML_DL {((1)) * ((K1()) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};
      IL_e_Bi[SiL]_WL -> IL_e_ML_ML {((1)) * ((K1()) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_e_Bi[SiL]_WL -> Bi[SL]_e_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(a1(u))) * (a1(v)))};
      IL_e_Bi[SiL]_WL -> IL_e_DL_DL {((1)) * ((K1()) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_e_Bi[SiL]_WL -> IL_e_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_e_Bi[SiL]_WL -> IL_e_DL_ML {((1)) * ((K1()) * ((1-b1(u))*(1-a1(u))) * (a1(v)))};
      IL_e_Bi[SiL]_WL -> Bi[SL]_e_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(1-a1(u))) * (1-a1(v)))};
      IL_e_Bi[SiL]_WL -> Bi[SL]_e_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(a1(u))) * (1-a1(v)))};

; IL_Bi[SiL]_e_WL (Bifurc)
;  (left):
      IL_Bi[SiL]_e_WL -> e_Si_e_e {(1)};
;  (center):
      IL_Bi[SiL]_e_WL -> IL_IL_e_WL {((1)) * ((b1(t)))};
      IL_Bi[SiL]_e_WL -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * (1))};
      IL_Bi[SiL]_e_WL -> Bi[SL]_Bp[eL]_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(v)))};
      IL_Bi[SiL]_e_WL -> Bi[SL]_Bp[eL]_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(v)))};
      IL_Bi[SiL]_e_WL -> IL_DL_e_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (a1(v)))};
      IL_Bi[SiL]_e_WL -> IL_ML_e_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (a1(v)))};
      IL_Bi[SiL]_e_WL -> Bi[SL]_Bm[SL]_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(v)))};
      IL_Bi[SiL]_e_WL -> IL_Bi[SiL]_e_WL {((1)) * ((b1(t)*pS()))};
      IL_Bi[SiL]_e_WL -> Bi[SL]_Bm[SL]_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(v)))};
      IL_Bi[SiL]_e_WL -> IL_ML_e_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (1-a1(v)))};
      IL_Bi[SiL]_e_WL -> IL_DL_e_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(v)))};

; IL_e_ML_Bi[SiL] (Bifurc)
;  (left):
      IL_e_ML_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_e_ML_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(u))*(1)) * ((1-b1(v))*(1)))};
      IL_e_ML_Bi[SiL] -> Bi[SL]_e_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_e_ML_Bi[SiL] -> IL_e_IL_WL {((1)) * ((b1(u)) * (1-b1(v)))};
      IL_e_ML_Bi[SiL] -> IL_e_ML_DL {((1)) * ((K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_e_ML_Bi[SiL] -> IL_e_ML_ML {((1)) * ((K1()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_e_ML_Bi[SiL] -> Bi[SL]_e_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_e_ML_Bi[SiL] -> IL_e_DL_DL {((1)) * ((K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_e_ML_Bi[SiL] -> IL_e_ML_Bi[SiL] {((1)) * ((b1(v)*pS()))};
      IL_e_ML_Bi[SiL] -> IL_e_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_e_ML_Bi[SiL] -> IL_e_DL_ML {((1)) * ((K1()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_e_ML_Bi[SiL] -> Bi[SL]_e_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_e_ML_Bi[SiL] -> Bi[SL]_e_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_e_ML_Bi[SiL] -> IL_e_ML_IL {((1)) * ((b1(v)))};

; IL_e_DL_Bi[SiL] (Bifurc)
;  (left):
      IL_e_DL_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_e_DL_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-g1(u))*(1)) * ((1-b1(v))*(1)))};
      IL_e_DL_Bi[SiL] -> Bi[SL]_e_Bp[eL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_e_DL_Bi[SiL] -> IL_e_DL_IL {((1)) * ((b1(v)))};
      IL_e_DL_Bi[SiL] -> IL_e_DL_Bi[SiL] {((1)) * ((b1(v)*pS()))};
      IL_e_DL_Bi[SiL] -> IL_e_IL_WL {((1)) * ((g1(u)) * (1-b1(v)))};
      IL_e_DL_Bi[SiL] -> IL_e_ML_DL {((1)) * ((K1()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_e_DL_Bi[SiL] -> IL_e_ML_ML {((1)) * ((K1()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_e_DL_Bi[SiL] -> Bi[SL]_e_Bm[SL]_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_e_DL_Bi[SiL] -> IL_e_DL_DL {((1)) * ((K1()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_e_DL_Bi[SiL] -> IL_e_Bi[SiL]_WL {((1)) * ((b1(u)*pS()))};
      IL_e_DL_Bi[SiL] -> IL_e_DL_ML {((1)) * ((K1()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(a1(v))))};
      IL_e_DL_Bi[SiL] -> Bi[SL]_e_Bp[eL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(u))*(1-a1(u))) * ((1-b1(v))*(1-a1(v))))};
      IL_e_DL_Bi[SiL] -> Bi[SL]_e_Bm[SL]_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(u))*(a1(u))) * ((1-b1(v))*(1-a1(v))))};

; IL_ML_e_Bi[SiL] (Bifurc)
;  (left):
      IL_ML_e_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_ML_e_Bi[SiL] -> IL_IL_e_WL {((1)) * ((b1(t)) * (1-b1(v)))};
      IL_ML_e_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * ((1-b1(v))*(1)))};
      IL_ML_e_Bi[SiL] -> IL_ML_e_IL {((1)) * ((b1(v)))};
      IL_ML_e_Bi[SiL] -> Bi[SL]_Bp[eL]_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v))))};
      IL_ML_e_Bi[SiL] -> Bi[SL]_Bp[eL]_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_e_Bi[SiL] -> IL_DL_e_ML {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v))))};
      IL_ML_e_Bi[SiL] -> IL_ML_e_ML {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v))))};
      IL_ML_e_Bi[SiL] -> Bi[SL]_Bm[SL]_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_e_Bi[SiL] -> IL_Bi[SiL]_e_WL {((1)) * ((b1(t)*pS()))};
      IL_ML_e_Bi[SiL] -> Bi[SL]_Bm[SL]_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(a1(v))))};
      IL_ML_e_Bi[SiL] -> IL_ML_e_DL {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v))))};
      IL_ML_e_Bi[SiL] -> IL_ML_e_Bi[SiL] {((1)) * ((b1(v)*pS()))};
      IL_ML_e_Bi[SiL] -> IL_DL_e_DL {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v))))};

; IL_ML_Bi[SiL]_e (Bifurc)
;  (left):
      IL_ML_Bi[SiL]_e -> e_e_Si_e {(1)};
;  (center):
      IL_ML_Bi[SiL]_e -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * ((1-b1(u))*(1)))};
      IL_ML_Bi[SiL]_e -> IL_ML_ML_e {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))))};
      IL_ML_Bi[SiL]_e -> IL_Bi[SiL]_WL_e {((1)) * ((b1(t)*pS()))};
      IL_ML_Bi[SiL]_e -> Bi[SL]_Bm[SL]_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(a1(u))))};
      IL_ML_Bi[SiL]_e -> IL_ML_IL_e {((1)) * ((b1(u)))};
      IL_ML_Bi[SiL]_e -> Bi[SL]_Bp[eL]_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))))};
      IL_ML_Bi[SiL]_e -> IL_IL_WL_e {((1)) * ((b1(t)) * (1-b1(u)))};
      IL_ML_Bi[SiL]_e -> Bi[SL]_Bp[eL]_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))))};
      IL_ML_Bi[SiL]_e -> IL_DL_DL_e {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))))};
      IL_ML_Bi[SiL]_e -> Bi[SL]_Bm[SL]_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))))};
      IL_ML_Bi[SiL]_e -> IL_ML_Bi[SiL]_e {((1)) * ((b1(u)*pS()))};
      IL_ML_Bi[SiL]_e -> IL_DL_ML_e {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))))};
      IL_ML_Bi[SiL]_e -> IL_ML_DL_e {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))))};

; IL_DL_e_Bi[SiL] (Bifurc)
;  (left):
      IL_DL_e_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_DL_e_Bi[SiL] -> IL_DL_e_IL {((1)) * ((b1(v)))};
      IL_DL_e_Bi[SiL] -> IL_IL_e_WL {((1)) * ((g1(t)) * (1-b1(v)))};
      IL_DL_e_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-g1(t))*(1)) * ((1-b1(v))*(1)))};
      IL_DL_e_Bi[SiL] -> Bi[SL]_Bp[eL]_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v))))};
      IL_DL_e_Bi[SiL] -> Bi[SL]_Bp[eL]_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_e_Bi[SiL] -> IL_DL_e_ML {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(a1(v))))};
      IL_DL_e_Bi[SiL] -> IL_ML_e_ML {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(a1(v))))};
      IL_DL_e_Bi[SiL] -> Bi[SL]_Bm[SL]_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_e_Bi[SiL] -> IL_Bi[SiL]_e_WL {((1)) * ((b1(t)*pS()))};
      IL_DL_e_Bi[SiL] -> Bi[SL]_Bm[SL]_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(a1(v))))};
      IL_DL_e_Bi[SiL] -> IL_DL_e_Bi[SiL] {((1)) * ((b1(v)*pS()))};
      IL_DL_e_Bi[SiL] -> IL_ML_e_DL {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(v))*(1-a1(v))))};
      IL_DL_e_Bi[SiL] -> IL_DL_e_DL {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(v))*(1-a1(v))))};

; IL_Bi[SiL]_WL_e (Bifurc)
;  (left):
      IL_Bi[SiL]_WL_e -> e_Si_e_e {(1)};
;  (center):
      IL_Bi[SiL]_WL_e -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)) * (1))};
      IL_Bi[SiL]_WL_e -> IL_ML_ML_e {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (a1(u)))};
      IL_Bi[SiL]_WL_e -> IL_Bi[SiL]_WL_e {((1)) * ((b1(t)*pS()))};
      IL_Bi[SiL]_WL_e -> Bi[SL]_Bm[SL]_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (a1(u)))};
      IL_Bi[SiL]_WL_e -> Bi[SL]_Bp[eL]_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (a1(u)))};
      IL_Bi[SiL]_WL_e -> IL_IL_WL_e {((1)) * ((b1(t)))};
      IL_Bi[SiL]_WL_e -> Bi[SL]_Bp[eL]_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)))};
      IL_Bi[SiL]_WL_e -> IL_DL_DL_e {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (1-a1(u)))};
      IL_Bi[SiL]_WL_e -> Bi[SL]_Bm[SL]_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))) * (1-a1(u)))};
      IL_Bi[SiL]_WL_e -> IL_DL_ML_e {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))) * (a1(u)))};
      IL_Bi[SiL]_WL_e -> IL_ML_DL_e {((1)) * ((K1()) * ((1-b1(t))*(a1(t))) * (1-a1(u)))};

; IL_DL_Bi[SiL]_e (Bifurc)
;  (left):
      IL_DL_Bi[SiL]_e -> e_e_Si_e {(1)};
;  (center):
      IL_DL_Bi[SiL]_e -> e_e_e_e {((1)) * ((1-K1()) * ((1-g1(t))*(1)) * ((1-b1(u))*(1)))};
      IL_DL_Bi[SiL]_e -> IL_ML_ML_e {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))))};
      IL_DL_Bi[SiL]_e -> IL_Bi[SiL]_WL_e {((1)) * ((b1(t)*pS()))};
      IL_DL_Bi[SiL]_e -> Bi[SL]_Bm[SL]_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(a1(u))))};
      IL_DL_Bi[SiL]_e -> IL_DL_IL_e {((1)) * ((b1(u)))};
      IL_DL_Bi[SiL]_e -> Bi[SL]_Bp[eL]_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))))};
      IL_DL_Bi[SiL]_e -> IL_IL_WL_e {((1)) * ((g1(t)) * (1-b1(u)))};
      IL_DL_Bi[SiL]_e -> Bi[SL]_Bp[eL]_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))))};
      IL_DL_Bi[SiL]_e -> IL_DL_DL_e {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(1-a1(u))))};
      IL_DL_Bi[SiL]_e -> Bi[SL]_Bm[SL]_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))))};
      IL_DL_Bi[SiL]_e -> IL_DL_ML_e {((1)) * ((K1()) * ((1-g1(t))*(1-a1(t))) * ((1-b1(u))*(a1(u))))};
      IL_DL_Bi[SiL]_e -> IL_DL_Bi[SiL]_e {((1)) * ((b1(u)*pS()))};
      IL_DL_Bi[SiL]_e -> IL_ML_DL_e {((1)) * ((K1()) * ((1-g1(t))*(a1(t))) * ((1-b1(u))*(1-a1(u))))};

; IL_e_e_Bi[SiL] (Bifurc)
;  (left):
      IL_e_e_Bi[SiL] -> e_e_e_Si {(1)};
;  (center):
      IL_e_e_Bi[SiL] -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(v))*(1)))};
      IL_e_e_Bi[SiL] -> IL_e_e_ML {((1)) * ((K1()) * ((1-b1(v))*(a1(v))))};
      IL_e_e_Bi[SiL] -> IL_e_e_DL {((1)) * ((K1()) * ((1-b1(v))*(1-a1(v))))};
      IL_e_e_Bi[SiL] -> Bi[SL]_e_e_Bm[SL] {((1)) * ((K1()*pS()) * ((1-b1(v))*(a1(v))))};
      IL_e_e_Bi[SiL] -> IL_e_e_IL {((1)) * ((b1(v)))};
      IL_e_e_Bi[SiL] -> Bi[SL]_e_e_Bp[eL] {((1)) * ((K1()*pS()) * ((1-b1(v))*(1-a1(v))))};
      IL_e_e_Bi[SiL] -> IL_e_e_Bi[SiL] {((1)) * ((b1(v)*pS()))};

; IL_e_Bi[SiL]_e (Bifurc)
;  (left):
      IL_e_Bi[SiL]_e -> e_e_Si_e {(1)};
;  (center):
      IL_e_Bi[SiL]_e -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(u))*(1)))};
      IL_e_Bi[SiL]_e -> IL_e_ML_e {((1)) * ((K1()) * ((1-b1(u))*(a1(u))))};
      IL_e_Bi[SiL]_e -> IL_e_Bi[SiL]_e {((1)) * ((b1(u)*pS()))};
      IL_e_Bi[SiL]_e -> Bi[SL]_e_Bp[eL]_e {((1)) * ((K1()*pS()) * ((1-b1(u))*(1-a1(u))))};
      IL_e_Bi[SiL]_e -> IL_e_IL_e {((1)) * ((b1(u)))};
      IL_e_Bi[SiL]_e -> Bi[SL]_e_Bm[SL]_e {((1)) * ((K1()*pS()) * ((1-b1(u))*(a1(u))))};
      IL_e_Bi[SiL]_e -> IL_e_DL_e {((1)) * ((K1()) * ((1-b1(u))*(1-a1(u))))};

; IL_Bi[SiL]_e_e (Bifurc)
;  (left):
      IL_Bi[SiL]_e_e -> e_Si_e_e {(1)};
;  (center):
      IL_Bi[SiL]_e_e -> e_e_e_e {((1)) * ((1-K1()) * ((1-b1(t))*(1)))};
      IL_Bi[SiL]_e_e -> Bi[SL]_Bp[eL]_e_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(1-a1(t))))};
      IL_Bi[SiL]_e_e -> Bi[SL]_Bm[SL]_e_e {((1)) * ((K1()*pS()) * ((1-b1(t))*(a1(t))))};
      IL_Bi[SiL]_e_e -> IL_DL_e_e {((1)) * ((K1()) * ((1-b1(t))*(1-a1(t))))};
      IL_Bi[SiL]_e_e -> IL_ML_e_e {((1)) * ((K1()) * ((1-b1(t))*(a1(t))))};
      IL_Bi[SiL]_e_e -> IL_IL_e_e {((1)) * ((b1(t)))};
      IL_Bi[SiL]_e_e -> IL_Bi[SiL]_e_e {((1)) * ((b1(t)*pS()))};

; Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] -> S_S_S_e {(1) * (1) * (1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_Bm[SL]_Bp[eL] -> L_L_L_L {(1) * (1) * (1) * (1)};

; Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] -> S_S_e_S {(1) * (1) * (1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_Bp[eL]_Bm[SL] -> L_L_L_L {(1) * (1) * (1) * (1)};

; Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] -> S_S_e_e {(1) * (1) * (1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_Bp[eL]_Bp[eL] -> L_L_L_L {(1) * (1) * (1) * (1)};

; Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] -> S_e_e_S {(1) * (1) * (1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_Bp[eL]_Bm[SL] -> L_L_L_L {(1) * (1) * (1) * (1)};

; Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] -> S_e_e_e {(1) * (1) * (1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_Bp[eL]_Bp[eL] -> L_L_L_L {(1) * (1) * (1) * (1)};

; Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] -> S_S_S_S {(1) * (1) * (1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_Bm[SL]_Bm[SL] -> L_L_L_L {(1) * (1) * (1) * (1)};

; Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] -> S_e_S_S {(1) * (1) * (1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_Bm[SL]_Bm[SL] -> L_L_L_L {(1) * (1) * (1) * (1)};

; Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] -> S_e_S_e {(1) * (1) * (1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_Bm[SL]_Bp[eL] -> L_L_L_L {(1) * (1) * (1) * (1)};

; Bi[SL]_e_Bm[SL]_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_e_Bm[SL]_Bp[eL] -> S_e_S_e {(1) * (1) * (1)};
;  (center):
      Bi[SL]_e_Bm[SL]_Bp[eL] -> L_e_L_L {(1) * (1) * (1)};

; Bi[SL]_e_Bp[eL]_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_e_Bp[eL]_Bp[eL] -> S_e_e_e {(1) * (1) * (1)};
;  (center):
      Bi[SL]_e_Bp[eL]_Bp[eL] -> L_e_L_L {(1) * (1) * (1)};

; Bi[SL]_e_Bm[SL]_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_e_Bm[SL]_Bm[SL] -> S_e_S_S {(1) * (1) * (1)};
;  (center):
      Bi[SL]_e_Bm[SL]_Bm[SL] -> L_e_L_L {(1) * (1) * (1)};

; Bi[SL]_e_Bp[eL]_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_e_Bp[eL]_Bm[SL] -> S_e_e_S {(1) * (1) * (1)};
;  (center):
      Bi[SL]_e_Bp[eL]_Bm[SL] -> L_e_L_L {(1) * (1) * (1)};

; Bi[SL]_Bp[eL]_e_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_e_Bm[SL] -> S_e_e_S {(1) * (1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_e_Bm[SL] -> L_L_e_L {(1) * (1) * (1)};

; Bi[SL]_Bm[SL]_e_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_e_Bp[eL] -> S_S_e_e {(1) * (1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_e_Bp[eL] -> L_L_e_L {(1) * (1) * (1)};

; Bi[SL]_Bm[SL]_e_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_e_Bm[SL] -> S_S_e_S {(1) * (1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_e_Bm[SL] -> L_L_e_L {(1) * (1) * (1)};

; Bi[SL]_Bp[eL]_e_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_e_Bp[eL] -> S_e_e_e {(1) * (1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_e_Bp[eL] -> L_L_e_L {(1) * (1) * (1)};

; Bi[SL]_Bp[eL]_Bp[eL]_e (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_Bp[eL]_e -> S_e_e_e {(1) * (1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_Bp[eL]_e -> L_L_L_e {(1) * (1) * (1)};

; Bi[SL]_Bm[SL]_Bp[eL]_e (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_Bp[eL]_e -> S_S_e_e {(1) * (1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_Bp[eL]_e -> L_L_L_e {(1) * (1) * (1)};

; Bi[SL]_Bm[SL]_Bm[SL]_e (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_Bm[SL]_e -> S_S_S_e {(1) * (1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_Bm[SL]_e -> L_L_L_e {(1) * (1) * (1)};

; Bi[SL]_Bp[eL]_Bm[SL]_e (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_Bm[SL]_e -> S_e_S_e {(1) * (1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_Bm[SL]_e -> L_L_L_e {(1) * (1) * (1)};

; Bi[SL]_e_e_Bm[SL] (Bifurc)
;  (left):
      Bi[SL]_e_e_Bm[SL] -> S_e_e_S {(1) * (1)};
;  (center):
      Bi[SL]_e_e_Bm[SL] -> L_e_e_L {(1) * (1)};

; Bi[SL]_e_e_Bp[eL] (Bifurc)
;  (left):
      Bi[SL]_e_e_Bp[eL] -> S_e_e_e {(1) * (1)};
;  (center):
      Bi[SL]_e_e_Bp[eL] -> L_e_e_L {(1) * (1)};

; Bi[SL]_e_Bp[eL]_e (Bifurc)
;  (left):
      Bi[SL]_e_Bp[eL]_e -> S_e_e_e {(1) * (1)};
;  (center):
      Bi[SL]_e_Bp[eL]_e -> L_e_L_e {(1) * (1)};

; Bi[SL]_e_Bm[SL]_e (Bifurc)
;  (left):
      Bi[SL]_e_Bm[SL]_e -> S_e_S_e {(1) * (1)};
;  (center):
      Bi[SL]_e_Bm[SL]_e -> L_e_L_e {(1) * (1)};

; Bi[SL]_Bp[eL]_e_e (Bifurc)
;  (left):
      Bi[SL]_Bp[eL]_e_e -> S_e_e_e {(1) * (1)};
;  (center):
      Bi[SL]_Bp[eL]_e_e -> L_L_e_e {(1) * (1)};

; Bi[SL]_Bm[SL]_e_e (Bifurc)
;  (left):
      Bi[SL]_Bm[SL]_e_e -> S_S_e_e {(1) * (1)};
;  (center):
      Bi[SL]_Bm[SL]_e_e -> L_L_e_e {(1) * (1)};

; Bi[SL]_e_e_e (Bifurc)
;  (left):
      Bi[SL]_e_e_e -> S_e_e_e {(1)};
;  (center):
      Bi[SL]_e_e_e -> L_e_e_e {(1)};

; e_e_e_Bi[SiLi] (Bifurc)
;  (left):
      e_e_e_Bi[SiLi] -> e_e_e_Si {(1)};
;  (center):
      e_e_e_Bi[SiLi] -> e_e_e_Li {(1)};

; e_e_Bi[SiLi]_e (Bifurc)
;  (left):
      e_e_Bi[SiLi]_e -> e_e_Si_e {(1)};
;  (center):
      e_e_Bi[SiLi]_e -> e_e_Li_e {(1)};

; e_Bi[SiLi]_e_e (Bifurc)
;  (left):
      e_Bi[SiLi]_e_e -> e_Si_e_e {(1)};
;  (center):
      e_Bi[SiLi]_e_e -> e_Li_e_e {(1)};


; Singlet machine:
;   7 total states
;   2 deterministic bifurcation states
; Branch machine:
;   21 total states
;   6 deterministic bifurcation states
; Composed machine:
;   287 = 225+61+1 effective total states
;   61 Bifurcation states (effectively deterministic)
;      (61 possibly non-deterministic bifurc states)
;   122 Emit states
;   81 Null states
;   225 (non-bifurcation) states in TM
;   686 transitions in TM
;      225^2 = 50625 (for comparison)
; Reduced composed machine:
;   230 = 168+61+1 effective total states
;   61 Bifurcation states (effectively deterministic)
;      (61 possibly non-deterministic bifurc states)
;   122 Emit states
;   24 Null states
;   168 (non-bifurcation) states in TM
;   1789 transitions in TM
;      168^2 = 28224 (for comparison)
; Fully-reduced composed machine:
;   540 = 144+395+1 effective total states
;   395 Bifurcation states (effectively deterministic)
;      (61 possibly non-deterministic bifurc states)
;   122 Emit states
;   0 Null states
;   144 (non-bifurcation) states in TM
;   1431 transitions in TM
;      144^2 = 20736 (for comparison)
