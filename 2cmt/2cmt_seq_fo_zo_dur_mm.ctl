;; ==========================================================================
;; 2CMT_SEQ_FO_ZO_DUR_MM
;;   Sequential absorption with first-order process with lag time followed by
;;     zero-order with duration
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --(zero-order, D2)--> [Central (V2)]
;;   (1-DF) x Dose at t=D2 --> [Depot] --Ka--^
;;   [Central (V2)] <--Q--> [Periph. (V3)]
;;   elimination from central: Vmax/Km
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
;;     CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at t=D2
;;                                  (ALAG1 = D2) then absorbed with Ka
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 2CMT seq_fo_zo_dur | Nonlinear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=3
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)
        COMP=(PERIPH)

; compartment 1 = depot
; compartment 2 = central (observation)
; compartment 3 = peripheral

$PK
; ---- typical values and between-subject variability ------------------
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V2    = THETA(3)*EXP(ETA(3))
  Q     = THETA(4)*EXP(ETA(4))
  V3    = THETA(5)*EXP(ETA(5))
  KA    = THETA(6)*EXP(ETA(6))
  D2    = THETA(7)*EXP(ETA(7))
  LGTDF = LOG(THETA(8)/(1-THETA(8))) + ETA(8)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)
  ET6   = ETA(6)
  ET7   = ETA(7)
  ET8   = ETA(8)

; ---- structural / input specification ---------------------------------
  F1    = 1 - DF
  F2    = DF
  ALAG1 = D2

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V2

$DES
  CONC    = A(2)/V2
  CP      = A(3)/V3
  ELR     = VM*CONC/(KM+CONC)
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - Q*CONC + Q*CP - ELR
  DADT(3) =  Q*CONC - Q*CP

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
$THETA  (0, 100.0)           ; 1 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 2 KM     Km (mg/L)
$THETA  (0, 50.0)            ; 3 V2     central volume (L)
$THETA  (0, 10.0)            ; 4 Q      inter-compartmental clearance (L/h)
$THETA  (0, 100.0)           ; 5 V3     peripheral volume (L)
$THETA  (0, 1.0)             ; 6 KA     absorption rate constant (1/h)
$THETA  (0, 2.0)             ; 7 D2     zero-order duration into cmt 2 (h)
$THETA  (0.001, 0.5, 0.999)  ; 8 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free it if the data support it)
$OMEGA  0.09             ; 3 IIV V2
$OMEGA  0.09             ; 4 IIV Q
$OMEGA  0.09             ; 5 IIV V3
$OMEGA  0.09             ; 6 IIV KA
$OMEGA  0 FIX            ; 7 IIV D2    (free it if the data support it)
$OMEGA  0 FIX            ; 8 IIV DF    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES PRED
         RES WRES ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8 SG1 SG2
         ONEHEADER NOPRINT FILE=2cmt_seq_fo_zo_dur_mm.tab
$TABLE   ID VM KM V2 Q V3 KA D2 DF ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8
         FIRSTONLY ONEHEADER NOPRINT FILE=2cmt_seq_fo_zo_dur_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=2cmt_seq_fo_zo_dur_mm_sim.tab

