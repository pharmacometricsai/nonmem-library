;; ==========================================================================
;; 3CMT_ZO_DUR_LIN
;; #13 of the three-compartment set
;;   Zero-order absorption with duration
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   Dose --(zero-order, D1)--> [Central (V1)]
;;   [Periph.1 (V2)] <--Q2--> [Central (V1)] <--Q3--> [Periph.2 (V3)]
;;   elimination from central: CL
;;
;; Data set / dosing requirements
;;   One dose record per administration: CMT=1 (central), RATE=-2
;;   (duration D1 modelled in $PK).
;;
;; Solver: closed-form (analytical)
;; ==========================================================================

$PROBLEM 3CMT zo_dur | Linear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN11 TRANS4

; compartment 1 = central (dosing + observation)
; compartment 2 = peripheral 1
; compartment 3 = peripheral 2
; compartment 4 = output

$PK
; ---- typical values and between-subject variability ------------------
  CL    = THETA(1)*EXP(ETA(1))
  V1    = THETA(2)*EXP(ETA(2))
  Q2    = THETA(3)*EXP(ETA(3))
  V2    = THETA(4)*EXP(ETA(4))
  Q3    = THETA(5)*EXP(ETA(5))
  V3    = THETA(6)*EXP(ETA(6))
  D1    = THETA(7)*EXP(ETA(7))

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
$THETA  (0, 10.0)            ; 3 Q2     inter-cmt clearance, periph.1 (L/h)
$THETA  (0, 100.0)           ; 4 V2     peripheral volume 1 (L)
$THETA  (0, 2.0)             ; 5 Q3     inter-cmt clearance, periph.2 (L/h)
$THETA  (0, 200.0)           ; 6 V3     peripheral volume 2 (L)
$THETA  (0, 2.0)             ; 7 D1     zero-order duration into cmt 1 (h)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V1
$OMEGA  0.09             ; 3 IIV Q2
$OMEGA  0.09             ; 4 IIV V2
$OMEGA  0.09             ; 5 IIV Q3
$OMEGA  0.09             ; 6 IIV V3
$OMEGA  0 FIX            ; 7 IIV D1    (free it if the data support it)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=3cmt_zo_dur_lin.tab
$TABLE   ID CL V1 Q2 V2 Q3 V3 D1 ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7
         FIRSTONLY ONEHEADER NOPRINT FILE=3cmt_zo_dur_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=3cmt_zo_dur_lin.tab

