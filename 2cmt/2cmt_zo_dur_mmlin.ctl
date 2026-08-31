;; ==========================================================================
;; 2CMT_ZO_DUR_MMLIN
;; #15 of the two-compartment set
;;   Zero-order absorption with duration
;;   Nonlinear and linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(zero-order, D1)--> [Central (V1)]
;;   [Central (V1)] <--Q--> [Periph. (V2)]
;;   elimination from central: CL + Vmax/Km
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (central), RATE=-2
;;   (duration D1 modelled in $PK).
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 2CMT zo_dur | Nonlinear and linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=2
        COMP=(CENTRAL, DEFDOSE, DEFOBS)
        COMP=(PERIPH)

; compartment 1 = central (dosing + observation)
; compartment 2 = peripheral

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V1    = THETA(2)*EXP(ETA(2))
  Q     = THETA(3)*EXP(ETA(3))
  V2    = THETA(4)*EXP(ETA(4))
  VM    = THETA(5)*EXP(ETA(5))
  KM    = THETA(6)*EXP(ETA(6))
  D1    = THETA(7)*EXP(ETA(7))

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V1

$DES
  CONC    = A(1)/V1
  CP      = A(2)/V2
  ELR     = CL*CONC + VM*CONC/(KM+CONC)
  DADT(1) = -Q*CONC + Q*CP - ELR
  DADT(2) =  Q*CONC - Q*CP

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     linear clearance (L/h)
$THETA  (0, 50.0)            ; 2 V1     central volume (L)
$THETA  (0, 10.0)            ; 3 Q      inter-compartmental clearance (L/h)
$THETA  (0, 100.0)           ; 4 V2     peripheral volume (L)
$THETA  (0, 100.0)           ; 5 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 6 KM     Km (mg/L)
$THETA  (0, 2.0)             ; 7 D1     zero-order duration into cmt 1 (h)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V1
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V2
$OMEGA  0.09             ; 5 IIV VM
$OMEGA  0 FIX            ; 6 IIV KM    (free it if the data support it)
$OMEGA  0 FIX            ; 7 IIV D1    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=2cmt_zo_dur_mmlin.tab
$TABLE   ID CL V1 Q V2 VM KM D1 ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7
         FIRSTONLY ONEHEADER NOPRINT FILE=2cmt_zo_dur_mmlin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=2cmt_zo_dur_mmlin.tab

