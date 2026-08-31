;; ==========================================================================
;; 2CMT_SIGMOID_DUR_LIN
;; #7 of the two-compartment set
;;   Sigmoid absorption with duration
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(zero-order, D1)--> [Depot] --Ka--> [Central (V2)]
;;   [Central (V2)] <--Q--> [Periph. (V3)]
;;   elimination from central: CL
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (depot), RATE=-2
;;   (duration D1 modelled in $PK).  Zero-order input INTO the depot
;;   followed by first-order transfer to central gives the sigmoid
;;   (S-shaped) plasma profile.
;;
;; Solver: closed-form (analytical)
;; ==========================================================================

$PROBLEM 2CMT sigmoid_dur | Linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN4 TRANS4

; compartment 1 = depot
; compartment 2 = central (observation)
; compartment 3 = peripheral
; compartment 4 = output

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V2    = THETA(2)*EXP(ETA(2))
  Q     = THETA(3)*EXP(ETA(3))
  V3    = THETA(4)*EXP(ETA(4))
  KA    = THETA(5)*EXP(ETA(5))
  D1    = THETA(6)*EXP(ETA(6))

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V2

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 5.0)             ; 1 CL     clearance (L/h)
$THETA  (0, 50.0)            ; 2 V2     central volume (L)
$THETA  (0, 10.0)            ; 3 Q      inter-compartmental clearance (L/h)
$THETA  (0, 100.0)           ; 4 V3     peripheral volume (L)
$THETA  (0, 1.0)             ; 5 KA     absorption rate constant (1/h)
$THETA  (0, 2.0)             ; 6 D1     zero-order duration into cmt 1 (h)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V2
$OMEGA  0.09             ; 3 IIV Q
$OMEGA  0.09             ; 4 IIV V3
$OMEGA  0.09             ; 5 IIV KA
$OMEGA  0 FIX            ; 6 IIV D1    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=2cmt_sigmoid_dur_lin.tab
$TABLE   ID CL V2 Q V3 KA D1 ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         FIRSTONLY ONEHEADER NOPRINT FILE=2cmt_sigmoid_dur_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=2cmt_sigmoid_dur_lin.tab

