;; ==========================================================================
;; 1CMT_SEQ_FO_ZO_DUR_MM
;; #29 of the one-compartment set
;;   Sequential absorption with first-order process with lag time followed by
;;     zero-order with duration
;;   Nonlinear elimination
;; --------------------------------------------------------------------------
;; Structure
;;   DF x Dose     --(zero-order, D2)--> [Central (V)] --Vmax/Km-->
;;   (1-DF) x Dose at t=D2 --> [Depot] --Ka--^
;;
;; Data set / dosing requirements
;;   TWO dose records per administration, both with the full AMT:
;;     CMT=2 (central), RATE=-2  -> DF     x Dose, zero-order over D2
;;     CMT=1 (depot),   RATE=0   -> (1-DF) x Dose, released at t=D2
;;                                  (ALAG1 = D2) then absorbed with Ka
;;
;; Solver: general non-linear ODE (ADVAN13)
;; ==========================================================================

$PROBLEM 1CMT seq_fo_zo_dur | Nonlinear elimination

$INPUT   ID TIME AMT RATE DV EVID MDV CMT ADDL II WT
$DATA    ../data/pkdata.csv IGNORE=@

$SUBROUTINES ADVAN13 TOL=9
$MODEL  NCOMPARTMENTS=2
        COMP=(DEPOT, DEFDOSE)
        COMP=(CENTRAL, DEFOBS)

; compartment 1 = depot
; compartment 2 = central (observation)

$PK
; ---- typical values and between-subject variability ------------------
  VM    = THETA(1)*EXP(ETA(1))
  KM    = THETA(2)*EXP(ETA(2))
  V     = THETA(3)*EXP(ETA(3))
  KA    = THETA(4)*EXP(ETA(4))
  D2    = THETA(5)*EXP(ETA(5))
  LGTDF = LOG(THETA(6)/(1-THETA(6))) + ETA(6)
  DF    = EXP(LGTDF)/(1+EXP(LGTDF))

; ---- structural / input specification ---------------------------------
  F1    = 1 - DF
  F2    = DF
  ALAG1 = D2

; ---- scaling (concentration = amount / volume) ------------------------
  S2 = V

$DES
  CONC    = A(2)/V
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - (VM*CONC/(KM+CONC))

$ERROR
  IPRED = F
  IF(IPRED.LE.0) IPRED = 1.0E-10
  Y     = IPRED + IPRED*EPS(1) + EPS(2)   ; combined prop. + add. error
  IRES  = DV - IPRED
  IWRES = IRES/IPRED

; ---- initial estimates -------------------------------------------------
$THETA  (0, 100.0)           ; 1 VM     Vmax (mg/h)
$THETA  (0, 5.0)             ; 2 KM     Km (mg/L)
$THETA  (0, 50.0)            ; 3 V      central volume (L)
$THETA  (0, 1.0)             ; 4 KA     absorption rate constant (1/h)
$THETA  (0, 2.0)             ; 5 D2     zero-order duration into cmt 2 (h)
$THETA  (0.001, 0.5, 0.999)  ; 6 DF     dose fraction (logit-scale IIV)

$OMEGA  0.09             ; 1 IIV VM
$OMEGA  0 FIX            ; 2 IIV KM    (free to estimate if supported)
$OMEGA  0.09             ; 3 IIV V
$OMEGA  0.09             ; 4 IIV KA
$OMEGA  0 FIX            ; 5 IIV D2    (free to estimate if supported)
$OMEGA  0 FIX            ; 6 IIV DF    (free to estimate if supported)

$SIGMA  0.04             ; 1 proportional residual error (CV ~20%)
$SIGMA  0.01             ; 2 additive residual error (variance)

$ESTIMATION METHOD=1 INTER MAXEVAL=9999 NSIG=3 SIGL=9 PRINT=5 NOABORT POSTHOC
$COVARIANCE PRINT=E UNCONDITIONAL

$TABLE   ID TIME AMT RATE EVID MDV CMT DV IPRED IRES IWRES CWRES
         PRED RES WRES
         ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_dur_mm.tab
$TABLE   ID VM KM V KA D2 DF ETA1 ETA2 ETA3 ETA4 ETA5 ETA6
         FIRSTONLY ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_dur_mm.par

;; ----------------------------------------------------------------------
;; SIMULATION BLOCK
;; To simulate instead of estimate: comment out $ESTIMATION and
;; $COVARIANCE above, then un-comment the two lines below.
;; ----------------------------------------------------------------------
; $SIMULATION (20260830) (20260831 NORMAL) ONLYSIM SUBPROBLEMS=200
; $TABLE ID TIME AMT RATE EVID MDV CMT DV IPRED
;        NOAPPEND ONEHEADER NOPRINT FILE=1cmt_seq_fo_zo_dur_mm.tab

