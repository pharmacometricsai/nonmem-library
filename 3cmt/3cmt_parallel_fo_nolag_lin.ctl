;; ==========================================================================
;; 3CMT_PARALLEL_FO_NOLAG_LIN
;; #22 of the three-compartment set
;;   Parallel first-order absorption without lag time
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --> [Depot1] --Ka1--\
;;                                      >--> [Central (V2)]
;;   (1-DF) x Dose --> [Depot2] --Ka2--/
;;   [Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
;;   elimination from central: CL
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=1 (depot 1), RATE=0   -> receives DF     x Dose via F1
;;     CMT=2 (depot 2), RATE=0   -> receives (1-DF) x Dose via F2
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 3CMT parallel_fo_nolag | Linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=5
        COMP=(DEPOT1, DEFDOSE)
        COMP=(DEPOT2)
        COMP=(CENTRAL, DEFOBS)
        COMP=(PERIPH1)
        COMP=(PERIPH2)

; compartment 1 = depot 1
; compartment 2 = depot 2
; compartment 3 = central (observation)
; compartment 4 = peripheral 1
; compartment 5 = peripheral 2

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q3    = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  Q4    = THETA(5)*EXP(ETA(5))
  V4    = THETA(6)*EXP(ETA(6))
  KA1   = THETA(7)*EXP(ETA(7))
  KA2   = THETA(8)*EXP(ETA(8))
  LGTDF = LOG(THETA(9)/(1-THETA(9))) + ETA(9)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- structural / input specification ---------------------------------
  F1 = DF
  F2 = 1 - DF

; ---- scaling (concentration = amount / volume) ------------------------
  S3 = V2

$DES
  CONC    = A(3)/V2
  CP1     = A(4)/V3
  CP2     = A(5)/V4
  DIST1   = Q3*(CONC - CP1)
  DIST2   = Q4*(CONC - CP2)
  ELR     = CL*CONC
  DADT(1) = -KA1*A(1)
  DADT(2) = -KA2*A(2)
  DADT(3) =  KA1*A(1) + KA2*A(2) - DIST1 - DIST2 - ELR
  DADT(4) =  DIST1
  DADT(5) =  DIST2

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     clearance (L/h)
$THETA  (0, 50.0)            ; 2 V2     central volume (L)
$THETA  (0, 10.0)            ; 3 Q3     inter-cmt clearance, periph.1 (L/h)
$THETA  (0, 100.0)           ; 4 V3     peripheral volume 1 (L)
$THETA  (0, 2.0)             ; 5 Q4     inter-cmt clearance, periph.2 (L/h)
$THETA  (0, 200.0)           ; 6 V4     peripheral volume 2 (L)
$THETA  (0, 1.5)             ; 7 KA1    absorption rate constant, depot 1 (1/h)
$THETA  (0, 0.3)             ; 8 KA2    absorption rate constant, depot 2 (1/h)
$THETA  (0.001, 0.5, 0.999)  ; 9 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q3
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV Q4
$OMEGA  0.09             ; 6 IIV V4
$OMEGA  0.09             ; 7 IIV KA1
$OMEGA  0.09             ; 8 IIV KA2
$OMEGA  0 FIX            ; 9 IIV DF    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=3cmt_parallel_fo_nolag_lin.tab
$TABLE   ID CL V2 Q3 V3 Q4 V4 KA1 KA2 DF ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         ETA7 ETA8 ETA9
         FIRSTONLY ONEHEADER NOPRINT FILE=3cmt_parallel_fo_nolag_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=3cmt_parallel_fo_nolag_lin.tab

