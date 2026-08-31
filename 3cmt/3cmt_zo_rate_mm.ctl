;; ==========================================================================
;; 3CMT_ZO_RATE_MM
;; #17 of the three-compartment set
;;   Zero-order absorption with rate
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(zero-order, R1)--> [Central (V1)]
;;   [Periph.1 (V2)] <--Q2--> [Central (V1)] <--Q3--> [Periph.2 (V3)]
;;   elimination from central: Vmax/Km
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (central), RATE=-1
;;   (rate R1 modelled in $PK).
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 3CMT zo_rate | Nonlinear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=3
        COMP=(CENTRAL, DEFDOSE, DEFOBS)
        COMP=(PERIPH1)
        COMP=(PERIPH2)

; compartment 1 = central (dosing + observation)
; compartment 2 = peripheral 1
; compartment 3 = peripheral 2

$PK
; ---- typical values and between-subject variability ------------------
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V1    = THETA(3)*EXP(ETA(3))
  Q2    = THETA(4)*EXP(ETA(4))
  V2    = THETA(5)*EXP(ETA(5))
  Q3    = THETA(6)*EXP(ETA(6))
  V3    = THETA(7)*EXP(ETA(7))
  R1    = THETA(8)*EXP(ETA(8))

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V1

$DES
  CONC    = A(1)/V1
  CP1     = A(2)/V2
  CP2     = A(3)/V3
  DIST1   = Q2*(CONC - CP1)
  DIST2   = Q3*(CONC - CP2)
  ELR     = VM*CONC/(KM+CONC)
  DADT(1) = -DIST1 - DIST2 - ELR
  DADT(2) =  DIST1
  DADT(3) =  DIST2

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 100.0)           ; 1 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 2 KM     Km (mg/L)
$THETA  (0, 50.0)            ; 3 V1     central volume (L)
$THETA  (0, 10.0)            ; 4 Q2     inter-cmt clearance, periph.1 (L/h)
$THETA  (0, 100.0)           ; 5 V2     peripheral volume 1 (L)
$THETA  (0, 2.0)             ; 6 Q3     inter-cmt clearance, periph.2 (L/h)
$THETA  (0, 200.0)           ; 7 V3     peripheral volume 2 (L)
$THETA  (0, 50.0)            ; 8 R1     zero-order rate into cmt 1 (mg/h)

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free it if the data support it)
$OMEGA  0.09             ; 3 IIV V1
$OMEGA  0.09             ; 4 IIV Q2
$OMEGA  0.09             ; 5 IIV V2
$OMEGA  0.09             ; 6 IIV Q3
$OMEGA  0.09             ; 7 IIV V3
$OMEGA  0 FIX            ; 8 IIV R1    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=3cmt_zo_rate_mm.tab
$TABLE   ID VM KM V1 Q2 V2 Q3 V3 R1 ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7
         ETA8
         FIRSTONLY ONEHEADER NOPRINT FILE=3cmt_zo_rate_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=3cmt_zo_rate_mm.tab

