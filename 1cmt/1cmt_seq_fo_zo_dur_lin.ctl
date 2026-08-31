;; ==========================================================================
;; 1CMT_SEQ_FO_ZO_DUR_LIN
;; #28 of the one-compartment set
;;   Sequential absorption with first-order process with lag time followed by
;;     zero-order with duration
;;   Linear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --(zero-order, D2)--> [Central (V)] --CL-->
;;   (1-DF) x Dose at t=D2 --> [Depot] --Ka--^
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
;;     CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at t=D2
;;                                  (ALAG1 = D2) then absorbed with Ka
;;
;; Solver: closed-form (analytical)
;; ==========================================================================

$PROBLEM 1CMT seq_fo_zo_dur | Linear elimination

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
  D2    = THETA(4)*EXP(ETA(4))
  LGTDF = LOG(THETA(5)/(1-THETA(5))) + ETA(5)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- structural / input specification ---------------------------------
  F1    = 1 - DF
  F2    = DF
  ALAG1 = D2

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
$THETA  (0, 2.0)             ; 4 D2     zero-order duration into cmt 2 (h)
$THETA  (0.001, 0.5, 0.999)  ; 5 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV CL
$OMEGA  0.09             ; 2 IIV V
$OMEGA  0.09             ; 3 IIV KA
$OMEGA  0 FIX            ; 4 IIV D2    (free to estimate if supported)
$OMEGA  0 FIX            ; 5 IIV DF    (free to estimate if supported)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_dur_lin.tab
$TABLE   ID CL V KA D2 DF ETA1 ETA2 ETA3 ETA4 ETA5
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_dur_lin.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_dur_lin.tab

