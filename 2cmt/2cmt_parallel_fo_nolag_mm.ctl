;; ==========================================================================
;; 2CMT_PARALLEL_FO_NOLAG_MM
;; #23 of the two-compartment set
;;   Parallel first-order absorption without lag time
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --> [Depot1] --Ka1--\
;;                                      >--> [Central (V2)]
;;   (1-DF) x Dose --> [Depot2] --Ka2--/
;;   [Central (V2)] <--Q--> [Periph. (V3)]
;;   elimination from central: Vmax/Km
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=1 (depot 1), RATE=0   -> receives DF     x Dose via F1
;;     CMT=2 (depot 2), RATE=0   -> receives (1-DF) x Dose via F2
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 2CMT parallel_fo_nolag | Nonlinear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=4
        COMP=(DEPOT1, DEFDOSE)
        COMP=(DEPOT2)
        COMP=(CENTRAL, DEFOBS)
        COMP=(PERIPH)

; compartment 1 = depot 1
; compartment 2 = depot 2
; compartment 3 = central (observation)
; compartment 4 = peripheral

$PK
; ---- typical values and between-subject variability ------------------
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V2    = THETA(3)*EXP(ETA(3))
  Q     = THETA(4)*EXP(ETA(4))
  V3    = THETA(5)*EXP(ETA(5))
  KA1   = THETA(6)*EXP(ETA(6))
  KA2   = THETA(7)*EXP(ETA(7))
  LGTDF = LOG(THETA(8)/(1-THETA(8))) + ETA(8)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- structural / input specification ---------------------------------
  F1 = DF
  F2 = 1 - DF

; ---- scaling (concentration = amount / volume) ------------------------
  S3 = V2

$DES
  CONC    = A(3)/V2
  CP      = A(4)/V3
  ELR     = VM*CONC/(KM+CONC)
  DADT(1) = -KA1*A(1)
  DADT(2) = -KA2*A(2)
  DADT(3) =  KA1*A(1) + KA2*A(2) - Q*CONC + Q*CP - ELR
  DADT(4) =  Q*CONC - Q*CP

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 100.0)           ; 1 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 2 KM     Km (mg/L)
$THETA  (0, 50.0)            ; 3 V2     central volume (L)
$THETA  (0, 10.0)            ; 4 Q      inter-compartmental clearance (L/h)
$THETA  (0, 100.0)           ; 5 V3     peripheral volume (L)
$THETA  (0, 1.5)             ; 6 KA1    absorption rate constant, depot 1 (1/h)
$THETA  (0, 0.3)             ; 7 KA2    absorption rate constant, depot 2 (1/h)
$THETA  (0.001, 0.5, 0.999)  ; 8 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free it if the data support it)
$OMEGA  0.09             ; 3 IIV V2
$OMEGA  0.09             ; 4 IIV Q
$OMEGA  0.09             ; 5 IIV V3
$OMEGA  0.09             ; 6 IIV KA1
$OMEGA  0.09             ; 7 IIV KA2
$OMEGA  0 FIX            ; 8 IIV DF    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=2cmt_parallel_fo_nolag_mm.tab
$TABLE   ID VM KM V2 Q V3 KA1 KA2 DF ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7
         ETA8
         FIRSTONLY ONEHEADER NOPRINT FILE=2cmt_parallel_fo_nolag_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=2cmt_parallel_fo_nolag_mm.tab

