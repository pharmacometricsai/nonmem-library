;; ==========================================================================
;; 3CMT_PARALLEL_FO_LAG_MM
;; #26 of the three-compartment set
;;   Parallel first-order absorption with lag time
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     at t=ALAG1 --> [Depot1] --Ka1--\
;;                                                 >--> [Central (V2)]
;;   (1-DF) x Dose at t=ALAG2 --> [Depot2] --Ka2--/
;;   [Periph.1 (V3)] <--Q3--> [Central (V2)] <--Q4--> [Periph.2 (V4)]
;;   elimination from central: Vmax/Km
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=1 (depot 1), RATE=0   -> DF     x Dose, delayed by ALAG1
;;     CMT=2 (depot 2), RATE=0   -> (1-DF) x Dose, delayed by ALAG2
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 3CMT parallel_fo_lag | Nonlinear elimination

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
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V2    = THETA(3)*EXP(ETA(3))
  Q3    = THETA(4)*EXP(ETA(4))
  V3    = THETA(5)*EXP(ETA(5))
  Q4    = THETA(6)*EXP(ETA(6))
  V4    = THETA(7)*EXP(ETA(7))
  KA1   = THETA(8)*EXP(ETA(8))
  KA2   = THETA(9)*EXP(ETA(9))
  LGTDF = LOG(THETA(10)/(1-THETA(10))) + ETA(10)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))
  ALAG1 = THETA(11)*EXP(ETA(11))
  ALAG2 = THETA(12)*EXP(ETA(12))

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
  ELR     = VM*CONC/(KM+CONC)
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
$THETA  (0, 100.0)           ; 1 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 2 KM     Km (mg/L)
$THETA  (0, 50.0)            ; 3 V2     central volume (L)
$THETA  (0, 10.0)            ; 4 Q3     inter-cmt clearance, periph.1 (L/h)
$THETA  (0, 100.0)           ; 5 V3     peripheral volume 1 (L)
$THETA  (0, 2.0)             ; 6 Q4     inter-cmt clearance, periph.2 (L/h)
$THETA  (0, 200.0)           ; 7 V4     peripheral volume 2 (L)
$THETA  (0, 1.5)             ; 8 KA1    absorption rate constant, depot 1 (1/h)
$THETA  (0, 0.3)             ; 9 KA2    absorption rate constant, depot 2 (1/h)
$THETA  (0.001, 0.5, 0.999)  ; 10 DF     dose fraction (logit-scale IIV)
$THETA  (0, 0.25)            ; 11 ALAG1  lag time, compartment 1 (h)
$THETA  (0, 1.0)             ; 12 ALAG2  lag time, compartment 2 (h)

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free it if the data support it)
$OMEGA  0.09             ; 3 IIV V2
$OMEGA  0.09             ; 4 IIV Q3
$OMEGA  0.09             ; 5 IIV V3
$OMEGA  0.09             ; 6 IIV Q4
$OMEGA  0.09             ; 7 IIV V4
$OMEGA  0.09             ; 8 IIV KA1
$OMEGA  0.09             ; 9 IIV KA2
$OMEGA  0 FIX            ; 10 IIV DF    (free it if the data support it)
$OMEGA  0 FIX            ; 11 IIV ALAG1 (free it if the data support it)
$OMEGA  0 FIX            ; 12 IIV ALAG2 (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=3cmt_parallel_fo_lag_mm.tab
$TABLE   ID VM KM V2 Q3 V3 Q4 V4 KA1 KA2 DF ALAG1 ALAG2 ETA1 ETA2 ETA3
         ETA4 ETA5 ETA6 ETA7 ETA8 ETA9 ETA10 ETA11 ETA12
         FIRSTONLY ONEHEADER NOPRINT FILE=3cmt_parallel_fo_lag_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=3cmt_parallel_fo_lag_mm.tab

