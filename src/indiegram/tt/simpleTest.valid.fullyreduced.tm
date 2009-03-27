;; Fully reduced transition matrix
; s_s_s_s (Start)
   s_s_s_s -> s_s_s_il  {(ilp(v))};
   s_s_s_s -> s_s_il_w  {(ilp(u)) * (1-ilp(v))};
   s_s_s_s -> s_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   s_s_s_s -> il_ml_ml_ml  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   s_s_s_s -> il_dl_ml_ml  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   s_s_s_s -> il_ml_dl_ml  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   s_s_s_s -> il_ml_ml_dl  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   s_s_s_s -> il_dl_dl_ml  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   s_s_s_s -> il_dl_ml_dl  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   s_s_s_s -> il_ml_dl_dl  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   s_s_s_s -> il_dl_dl_dl  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   s_s_s_s -> e_e_e_e  {(ep()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; s_s_s_il (Emit: 4)
   s_s_s_il -> s_s_s_il  {(ilp(v))};
   s_s_s_il -> s_s_il_w  {(ilp(u)) * (1-ilp(v))};
   s_s_s_il -> s_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   s_s_s_il -> il_ml_ml_ml  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   s_s_s_il -> il_dl_ml_ml  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   s_s_s_il -> il_ml_dl_ml  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   s_s_s_il -> il_ml_ml_dl  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   s_s_s_il -> il_dl_dl_ml  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   s_s_s_il -> il_dl_ml_dl  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   s_s_s_il -> il_ml_dl_dl  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   s_s_s_il -> il_dl_dl_dl  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   s_s_s_il -> e_e_e_e  {(ep()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; s_s_il_w (Emit: 3)
   s_s_il_w -> s_s_il_w  {(ilp(u))};
   s_s_il_w -> s_il_w_w  {(ilp(t)) * (1-ilp(u))};
   s_s_il_w -> il_ml_ml_ml  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * (mlp(v))};
   s_s_il_w -> il_dl_ml_ml  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * (mlp(v))};
   s_s_il_w -> il_ml_dl_ml  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (mlp(v))};
   s_s_il_w -> il_ml_ml_dl  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * (1-mlp(v))};
   s_s_il_w -> il_dl_dl_ml  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (mlp(v))};
   s_s_il_w -> il_dl_ml_dl  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * (1-mlp(v))};
   s_s_il_w -> il_ml_dl_dl  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (1-mlp(v))};
   s_s_il_w -> il_dl_dl_dl  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (1-mlp(v))};
   s_s_il_w -> e_e_e_e  {(ep()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * (1)};
; s_il_w_w (Emit: 2)
   s_il_w_w -> s_il_w_w  {(ilp(t))};
   s_il_w_w -> il_ml_ml_ml  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * (mlp(u)) * (mlp(v))};
   s_il_w_w -> il_dl_ml_ml  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * (mlp(u)) * (mlp(v))};
   s_il_w_w -> il_ml_dl_ml  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * (1-mlp(u)) * (mlp(v))};
   s_il_w_w -> il_ml_ml_dl  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * (mlp(u)) * (1-mlp(v))};
   s_il_w_w -> il_dl_dl_ml  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * (1-mlp(u)) * (mlp(v))};
   s_il_w_w -> il_dl_ml_dl  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * (mlp(u)) * (1-mlp(v))};
   s_il_w_w -> il_ml_dl_dl  {(1-ep()) * ((1-ilp(t))*(mlp(t))) * (1-mlp(u)) * (1-mlp(v))};
   s_il_w_w -> il_dl_dl_dl  {(1-ep()) * ((1-ilp(t))*(1-mlp(t))) * (1-mlp(u)) * (1-mlp(v))};
   s_il_w_w -> e_e_e_e  {(ep()) * ((1-ilp(t))*(1)) * (1) * (1)};
; il_ml_ml_il (Emit: 4)
   il_ml_ml_il -> il_ml_ml_il  {(ilp(v))};
   il_ml_ml_il -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_il -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_il -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_il -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_il -> il_ml_il_w  {(ilp(u)) * (1-ilp(v))};
   il_ml_ml_il -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_il -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_il -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_il -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_il -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_ml_ml_il -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_ml_ml_il -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_dl_ml_il (Emit: 4)
   il_dl_ml_il -> il_dl_ml_il  {(ilp(v))};
   il_dl_ml_il -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_il -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_il -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_il -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_il -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_il -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_il -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_il -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_il -> il_dl_il_w  {(ilp(u)) * (1-ilp(v))};
   il_dl_ml_il -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_dl_ml_il -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_dl_ml_il -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_ml_dl_il (Emit: 4)
   il_ml_dl_il -> il_ml_dl_il  {(ilp(v))};
   il_ml_dl_il -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_il -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_il -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_il -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_il -> il_ml_il_w  {(ilp(u)) * (1-ilp(v))};
   il_ml_dl_il -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_il -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_il -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_il -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_il -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_ml_dl_il -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_ml_dl_il -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_ml_ml_ml (Emit: 1,2,3,4)
   il_ml_ml_ml -> il_ml_ml_il  {(ilp(v))};
   il_ml_ml_ml -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_ml -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_ml -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_ml -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_ml -> il_ml_il_w  {(ilp(u)) * (1-ilp(v))};
   il_ml_ml_ml -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_ml -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_ml -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_ml -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_ml -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_ml_ml_ml -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_ml_ml_ml -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_dl_dl_il (Emit: 4)
   il_dl_dl_il -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_il -> il_dl_dl_il  {(ilp(v))};
   il_dl_dl_il -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_il -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_il -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_il -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_il -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_il -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_il -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_il -> il_dl_il_w  {(ilp(u)) * (1-ilp(v))};
   il_dl_dl_il -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_dl_dl_il -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_dl_dl_il -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_dl_ml_ml (Emit: 1,2,3,4)
   il_dl_ml_ml -> il_dl_ml_il  {(ilp(v))};
   il_dl_ml_ml -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_ml -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_ml -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_ml -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_ml -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_ml -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_ml -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_ml -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_ml -> il_dl_il_w  {(ilp(u)) * (1-ilp(v))};
   il_dl_ml_ml -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_dl_ml_ml -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_dl_ml_ml -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_ml_dl_ml (Emit: 1,2,3,4)
   il_ml_dl_ml -> il_ml_dl_il  {(ilp(v))};
   il_ml_dl_ml -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_ml -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_ml -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_ml -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_ml -> il_ml_il_w  {(ilp(u)) * (1-ilp(v))};
   il_ml_dl_ml -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_ml -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_ml -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_ml -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_ml -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_ml_dl_ml -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_ml_dl_ml -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_ml_ml_dl (Emit: 1,2,3,4)
   il_ml_ml_dl -> il_ml_ml_il  {(ilp(v))};
   il_ml_ml_dl -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_dl -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_dl -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_dl -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_dl -> il_ml_il_w  {(ilp(u)) * (1-ilp(v))};
   il_ml_ml_dl -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_ml_dl -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_dl -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_dl -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_ml_dl -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_ml_ml_dl -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_ml_ml_dl -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_dl_dl_ml (Emit: 1,2,3,4)
   il_dl_dl_ml -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_ml -> il_dl_dl_il  {(ilp(v))};
   il_dl_dl_ml -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_ml -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_ml -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_ml -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_ml -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_ml -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_ml -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_ml -> il_dl_il_w  {(ilp(u)) * (1-ilp(v))};
   il_dl_dl_ml -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_dl_dl_ml -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_dl_dl_ml -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_dl_ml_dl (Emit: 1,2,3,4)
   il_dl_ml_dl -> il_dl_ml_il  {(ilp(v))};
   il_dl_ml_dl -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_dl -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_dl -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_dl -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_dl -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_ml_dl -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_dl -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_dl -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_ml_dl -> il_dl_il_w  {(ilp(u)) * (1-ilp(v))};
   il_dl_ml_dl -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_dl_ml_dl -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_dl_ml_dl -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_ml_dl_dl (Emit: 1,2,3,4)
   il_ml_dl_dl -> il_ml_dl_il  {(ilp(v))};
   il_ml_dl_dl -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_dl -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_dl -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_dl -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_dl -> il_ml_il_w  {(ilp(u)) * (1-ilp(v))};
   il_ml_dl_dl -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_ml_dl_dl -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_dl -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_dl -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_ml_dl_dl -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_ml_dl_dl -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_ml_dl_dl -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_ml_il_w (Emit: 3)
   il_ml_il_w -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * (mlp(v))};
   il_ml_il_w -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * (mlp(v))};
   il_ml_il_w -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (mlp(v))};
   il_ml_il_w -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * (1-mlp(v))};
   il_ml_il_w -> il_ml_il_w  {(ilp(u))};
   il_ml_il_w -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (mlp(v))};
   il_ml_il_w -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * (1-mlp(v))};
   il_ml_il_w -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (1-mlp(v))};
   il_ml_il_w -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (1-mlp(v))};
   il_ml_il_w -> il_il_w_w  {(ilp(t)) * (1-ilp(u))};
   il_ml_il_w -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * (1)};
   il_ml_il_w -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * (1)};
; il_dl_dl_dl (Emit: 1,2,3,4)
   il_dl_dl_dl -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_dl -> il_dl_dl_il  {(ilp(v))};
   il_dl_dl_dl -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_dl -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_dl -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_dl -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(mlp(v)))};
   il_dl_dl_dl -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_dl -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_dl -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * ((1-ilp(v))*(1-mlp(v)))};
   il_dl_dl_dl -> il_dl_il_w  {(ilp(u)) * (1-ilp(v))};
   il_dl_dl_dl -> il_il_w_w  {(ilp(t)) * (1-ilp(u)) * (1-ilp(v))};
   il_dl_dl_dl -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
   il_dl_dl_dl -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * ((1-ilp(v))*(1))};
; il_dl_il_w (Emit: 3)
   il_dl_il_w -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * (mlp(v))};
   il_dl_il_w -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * (mlp(v))};
   il_dl_il_w -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (mlp(v))};
   il_dl_il_w -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(mlp(u))) * (1-mlp(v))};
   il_dl_il_w -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (mlp(v))};
   il_dl_il_w -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(mlp(u))) * (1-mlp(v))};
   il_dl_il_w -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (1-mlp(v))};
   il_dl_il_w -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * ((1-ilp(u))*(1-mlp(u))) * (1-mlp(v))};
   il_dl_il_w -> il_dl_il_w  {(ilp(u))};
   il_dl_il_w -> il_il_w_w  {(ilp(t)) * (1-ilp(u))};
   il_dl_il_w -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * (1)};
   il_dl_il_w -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * ((1-ilp(u))*(1)) * (1)};
; il_il_w_w (Emit: 2)
   il_il_w_w -> il_ml_ml_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * (mlp(u)) * (mlp(v))};
   il_il_w_w -> il_dl_ml_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * (mlp(u)) * (mlp(v))};
   il_il_w_w -> il_ml_dl_ml  {(ilp()) * ((1-ilp(t))*(mlp(t))) * (1-mlp(u)) * (mlp(v))};
   il_il_w_w -> il_ml_ml_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * (mlp(u)) * (1-mlp(v))};
   il_il_w_w -> il_dl_dl_ml  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * (1-mlp(u)) * (mlp(v))};
   il_il_w_w -> il_dl_ml_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * (mlp(u)) * (1-mlp(v))};
   il_il_w_w -> il_ml_dl_dl  {(ilp()) * ((1-ilp(t))*(mlp(t))) * (1-mlp(u)) * (1-mlp(v))};
   il_il_w_w -> il_dl_dl_dl  {(ilp()) * ((1-ilp(t))*(1-mlp(t))) * (1-mlp(u)) * (1-mlp(v))};
   il_il_w_w -> il_il_w_w  {(ilp(t))};
   il_il_w_w -> Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss]  {(blp()) * ((1-ilp(t))*(1)) * (1) * (1)};
   il_il_w_w -> e_e_e_e  {(1-ilp()-blp()) * ((1-ilp(t))*(1)) * (1) * (1)};

;; Bifurcations
; Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss] (Bifurc)
;  (left):
      Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss] -> s_s_s_s {(1) * (1) * (1) * (1)};
;  (center):
      Bi[ss]_Bm[ss]_Bm[ss]_Bm[ss] -> s_s_s_s {(1) * (1) * (1) * (1)};

; Some statistics on the composed machine(s)

; Singlet machine:
;   4 total states
;   1 deterministic bifurcation states
; Branch machine:
;   7 total states
;   1 deterministic bifurcation states
; Composed machine:
;   31 = 29+1+1 total states
;      4*7*7*7 = 1372 (for comparison)
;   1 deterministic bifurc states
;   29 (non-bifurcation) states in TM
;   73 transitions in TM
;      29^2 = 841 (for comparison)
; Reduced composed machine:
;   21 = 19+1+1 total states
;   1 deterministic bifurc states
;   19 (non-bifurcation) states in reduced TM
;   236 transitions in TM
;      19^2 = 361 (for comparison)
; Fully reduced composed machine:
;   21 = 19+1+1 effective total states
;   1 effective deterministic bifurcation states
;    (1 possibly non-deterministic bifurc states)
;   19 (non-bifurcation) states in TM
;   236 transitions in reduced TM
;      19^2 = 361 (for comparison)
