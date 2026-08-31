;; ==========================================================================
;; 2CMT_SEQ_FO_ZO_DUR_LIN
;;   Sequential absorption with first-order process with lag time followed by
;;     zero-order with duration
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --(zero-order, D2)--> [Central (V2)]
;;   (1-DF) x Dose at t=D2 --> [Depot] --Ka--^
;;   [Central (V2)] <--Q--> [Periph. (V3)]
;;   elimination from central: CL
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
;;     CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at t=D2
;;                                  (ALAG1 = D2) then absorbed with Ka
;;
;; Solver: closed-form (analytical)
;; ==========================================================================

$PROBLEM 2CMT seq_fo_zo_dur | Linear elimination

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
  D2    = THETA(6)*EXP(ETA(6))
  LGTDF = LOG(THETA(7)/(1-THETA(7))) + ETA(7)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)
  ET6   = ETA(6)
  ET7   = ETA(7)

; ---- structural / input specification ---------------------------------
  F1    = 1 - DF
  F2    = DF
  ALAG1 = D2

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V2

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- residual error variances exported to $TABLE ----------------------
  SG1   = SIGMA(1,1)          ; proportional
  SG2   = SIGMA(2,2)          ; additive

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     clearance (L/h)
$THETA  (0, 50.0)            ; 2 V2     central volume (L)
$THETA  (0, 10.0)            ; 3 Q      inter-compartmental clearance (L/h)
$THETA  (0, 100.0)           ; 4 V3     peripheral volume (L)
$THETA  (0, 1.0)             ; 5 KA     absorption rate constant (1/h)
$THETA  (0, 2.0)             ; 6 D2     zero-order duration into cmt 2 (h)
$THETA  (0.001, 0.5, 0.999)  ; 7 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV KA
$OMEGA  0 FIX            ; 6 IIV D2    (free it if the data support it)
$OMEGA  0 FIX            ; 7 IIV DF    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES PRED
         RES WRES ET1 ET2 ET3 ET4 ET5 ET6 ET7 SG1 SG2
         ONEHEADER NOPRINT FILE=2cmt_seq_fo_zo_dur_lin.tab
$TABLE   ID CL V2 Q V3 KA D2 DF ET1 ET2 ET3 ET4 ET5 ET6 ET7
         FIRSTONLY ONEHEADER NOPRINT FILE=2cmt_seq_fo_zo_dur_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=2cmt_seq_fo_zo_dur_lin_sim.tab

