;; ==========================================================================
;; 2CMT_SIGMOID_RATE_MMLIN
;;   Sigmoid absorption with rate
;;   Nonlinear and linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(zero-order, R1)--> [Depot] --Ka--> [Central (V2)]
;;   [Central (V2)] <--Q--> [Periph. (V3)]
;;   elimination from central: CL + Vmax/Km
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (depot), RATE=-1
;;   (rate R1 modelled in $PK).
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 2CMT sigmoid_rate | Nonlinear and linear elimination

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
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q     = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  VM    = THETA(5)*EXP(ETA(5))
  KM    = THETA(6)*EXP(ETA(6))
  KA    = THETA(7)*EXP(ETA(7))
  R1    = THETA(8)*EXP(ETA(8))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)
  ET5   = ETA(5)
  ET6   = ETA(6)
  ET7   = ETA(7)
  ET8   = ETA(8)

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V2

$DES
  CONC    = A(2)/V2
  CP      = A(3)/V3
  ELR     = CL*CONC + VM*CONC/(KM+CONC)
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
$THETA  (0, 5.0)             ; 1 CL     linear clearance (L/h)
$THETA  (0, 50.0)            ; 2 V2     central volume (L)
$THETA  (0, 10.0)            ; 3 Q      inter-compartmental clearance (L/h)
$THETA  (0, 100.0)           ; 4 V3     peripheral volume (L)
$THETA  (0, 100.0)           ; 5 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 6 KM     Km (mg/L)
$THETA  (0, 1.0)             ; 7 KA     absorption rate constant (1/h)
$THETA  (0, 50.0)            ; 8 R1     zero-order rate into cmt 1 (mg/h)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV VM
$OMEGA  0 FIX            ; 6 IIV KM    (free it if the data support it)
$OMEGA  0.09             ; 7 IIV KA
$OMEGA  0 FIX            ; 8 IIV R1    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES PRED
         RES WRES ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8 SG1 SG2
         ONEHEADER NOPRINT FILE=2cmt_sigmoid_rate_mmlin.tab
$TABLE   ID CL V2 Q V3 VM KM KA R1 ET1 ET2 ET3 ET4 ET5 ET6 ET7 ET8
         FIRSTONLY ONEHEADER NOPRINT FILE=2cmt_sigmoid_rate_mmlin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=2cmt_sigmoid_rate_mmlin_sim.tab

