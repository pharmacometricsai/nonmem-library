;; ==========================================================================
;; 1CMT_ZO_DUR_MM
;;   Zero-order absorption with duration
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(zero-order, D1)--> [Central (V)] --Vmax/Km-->
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (central), RATE=-2
;;   (duration D1 modelled in $PK).
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 1CMT zo_dur | Nonlinear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=1
        COMP=(CENTRAL, DEFDOSE, DEFOBS)

; compartment 1 = central (dosing + observation)

$PK
; ---- typical values and between-subject variability ------------------
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V     = THETA(3)*EXP(ETA(3))
  D1    = THETA(4)*EXP(ETA(4))

; ---- random effects exported to $TABLE --------------------------------
  ET1   = ETA(1)
  ET2   = ETA(2)
  ET3   = ETA(3)
  ET4   = ETA(4)

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V

$DES
  CONC    = A(1)/V
  DADT(1) = -(VM*CONC/(KM+CONC))

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
$THETA  (0, 50.0)            ; 3 V      central volume (L)
$THETA  (0, 2.0)             ; 4 D1     zero-order duration into cmt 1 (h)

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free to estimate if supported)
$OMEGA  0.09             ; 3 IIV V
$OMEGA  0 FIX            ; 4 IIV D1    (free to estimate if supported)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES PRED
         RES WRES ET1 ET2 ET3 ET4 SG1 SG2
         ONEHEADER NOPRINT FILE=1cmt_zo_dur_mm.tab
$TABLE   ID VM KM V D1 ET1 ET2 ET3 ET4
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_zo_dur_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_zo_dur_mm_sim.tab

