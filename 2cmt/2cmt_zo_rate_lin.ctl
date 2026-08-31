;; ==========================================================================
;; 2CMT_ZO_RATE_LIN
;; #16 of the two-compartment set
;;   Zero-order absorption with rate
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(zero-order, R1)--> [Central (V1)]
;;   [Central (V1)] <--Q--> [Periph. (V2)]
;;   elimination from central: CL
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (central), RATE=-1
;;   (rate R1 modelled in $PK).
;;
;; Solver: closed-form (analytical)
;; ==========================================================================

$PROBLEM 2CMT zo_rate | Linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN3 TRANS4

; compartment 1 = central (dosing + observation)
; compartment 2 = peripheral
; compartment 3 = output

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V1    = THETA(2)*EXP(ETA(2))
  Q     = THETA(3)*EXP(ETA(3))
  V2    = THETA(4)*EXP(ETA(4))
  R1    = THETA(5)*EXP(ETA(5))

; ---- scaling (concentration = amount / volume) ------------------------
  S1 = V1

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     clearance (L/h)
$THETA  (0, 50.0)            ; 2 V1     central volume (L)
$THETA  (0, 10.0)            ; 3 Q      inter-compartmental clearance (L/h)
$THETA  (0, 100.0)           ; 4 V2     peripheral volume (L)
$THETA  (0, 50.0)            ; 5 R1     zero-order rate into cmt 1 (mg/h)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V1
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V2
$OMEGA  0 FIX            ; 5 IIV R1    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=2cmt_zo_rate_lin.tab
$TABLE   ID CL V1 Q V2 R1 ETA1 ETA2 ETA3 ETA4 ETA5
         FIRSTONLY ONEHEADER NOPRINT FILE=2cmt_zo_rate_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=2cmt_zo_rate_lin.tab

