;; ==========================================================================
;; 1CMT_FO_NOLAG_LIN
;; #1 of the one-compartment set
;;   First-order absorption without lag time
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(instantaneous)--> [Depot] --Ka--> [Central (V)] --CL-->
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (depot), RATE=0.
;;
;; Solver: closed-form (analytical)
;; ==========================================================================

$PROBLEM 1CMT fo_nolag | Linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN2 TRANS2

; compartment 1 = depot
; compartment 2 = central (observation)
; compartment 3 = output

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V     = THETA(2)*EXP(ETA(2))
  KA    = THETA(3)*EXP(ETA(3))

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     clearance (L/h)
$THETA  (0, 50.0)            ; 2 V      central volume (L)
$THETA  (0, 1.0)             ; 3 KA     absorption rate constant (1/h)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV KA

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=1cmt_fo_nolag_lin.tab
$TABLE   ID CL V KA ETA1 ETA2 ETA3
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_fo_nolag_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_fo_nolag_lin.tab

