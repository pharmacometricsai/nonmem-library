;; ==========================================================================
;; 2CMT_PARALLEL_FO_ZO_RATE_LIN
;; #37 of the two-compartment set
;;   Parallel absorption with first-order process and zero-order with rate
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --(zero-order, R2)--> [Central (V2)]
;;   (1-DF) x Dose --> [Depot] --Ka--^
;;   [Central (V2)] <--Q--> [Periph. (V3)]
;;   elimination from central: CL
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=2 (central), RATE=-1  -> DF     x Dose, zero-order at rate R2
;;     CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, absorbed with Ka
;;
;; Solver: closed-form (analytical)
;; ==========================================================================

$PROBLEM 2CMT parallel_fo_zo_rate | Linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN4 TRANS4

; compartment 1 = depot
; compartment 2 = central (observation)
; compartment 3 = peripheral
; compartment 4 = output

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q     = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  KA    = THETA(5)*EXP(ETA(5))
  R2    = THETA(6)*EXP(ETA(6))
  LGTDF = LOG(THETA(7)/(1-THETA(7))) + ETA(7)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- structural / input specification ---------------------------------
  F1 = 1 - DF
  F2 = DF

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V2

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     clearance (L/h)
$THETA  (0, 50.0)            ; 2 V2     central volume (L)
$THETA  (0, 10.0)            ; 3 Q      inter-compartmental clearance (L/h)
$THETA  (0, 100.0)           ; 4 V3     peripheral volume (L)
$THETA  (0, 1.0)             ; 5 KA     absorption rate constant (1/h)
$THETA  (0, 50.0)            ; 6 R2     zero-order rate into cmt 2 (mg/h)
$THETA  (0.001, 0.5, 0.999)  ; 7 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV KA
$OMEGA  0 FIX            ; 6 IIV R2    (free it if the data support it)
$OMEGA  0 FIX            ; 7 IIV DF    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=2cmt_parallel_fo_zo_rate_lin.tab
$TABLE   ID CL V2 Q V3 KA R2 DF ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7
         FIRSTONLY ONEHEADER NOPRINT FILE=2cmt_parallel_fo_zo_rate_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=2cmt_parallel_fo_zo_rate_lin.tab

